      subroutine gwflow_pump_allo(ob_id_dmd,ob_id_src,demand_vol,extracted,dmd_unmet) !rtb gwflow

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine removes groundwater from gwflow grid cells to meet an irrigation demand
!!    (called from actions (irrigate from "gwf") and wallo_withdraw (POD type "gwf");
!!    the pumped volumes are used in gwflow_simulate, in the groundwater balance equations)
!!
!!    source cells:
!!      ob_id_src = 0  the cells intersected by HRU ob_id_dmd, with the demand split by intersected area
!!      ob_id_src > 0  that one cell (cell id as numbered in the gwflow input)
!!    each cell supplies at most its unconfined storage, (head - bottom) * area * Sy, less the pumping
!!    already taken from it on the same day; demand left unmet after the first split is taken from the
!!    cells that still have storage, in proportion to that storage
!!    delivered and unmet volumes are summed per HRU (hru_pump, hru_pump_unmet, demand in hru_pump_dmd)
!!    and per cell (gw_hyd_ss%ppag, gw_hyd_ss%ppdf)

      use gwflow_module
      use organic_mineral_mass_module, only : soil1
      use hru_module, only : hru
      use hydrograph_module, only : wet
      use constituent_mass_module
      use res_salt_module, only : wetsalt_d
      use res_cs_module, only : wetcs_d
      use salt_module, only : hsaltb_d
      use cs_module, only : hcsb_d

      implicit none

      integer, intent (in) :: ob_id_dmd     !       |HRU with the demand (0 = no HRU, e.g. water allocation without HRUs)
      integer, intent (in) :: ob_id_src     !       |source cell (0 = cells intersected by the HRU)
      real, intent (in) :: demand_vol       !m3     |volume of water demand
      real, intent (inout) :: extracted     !m3     |volume of groundwater extracted from the aquifer
      real, intent (inout) :: dmd_unmet     !m3     |volume of demand not met by the aquifer
      integer :: i                          !       |counter
      integer :: ncl                        !       |number of source cells
      integer :: cell_id                    !       |gwflow cell
      integer, dimension (:), allocatable :: cells  !    |source cells
      real, dimension (:), allocatable :: wts       !    |fraction of the demand assigned to each cell
      real, dimension (:), allocatable :: avail     !m3  |groundwater available in each cell
      real, dimension (:), allocatable :: req       !m3  |demand assigned to each cell
      real, dimension (:), allocatable :: removed   !m3  |groundwater removed from each cell
      real :: sum_wts                       !       |sum of the intersected area fractions
      real :: remaining                     !m3     |demand not met after the first split
      real :: spare                         !m3     |storage left in all cells after the first split
      real :: add                           !m3     |volume added to a cell
      real :: shortfall                     !m3     |sum of cell shortfalls (assigned demand - removed)

      extracted = 0.
      dmd_unmet = 0.
      if (demand_vol <= 0.) return

      !source cells and the fraction of the demand assigned to each
      ncl = 0
      if (ob_id_src > 0) then
        ncl = 1
        allocate (cells(1), wts(1))
        if (grid_type == "structured") then
          cells(1) = cell_id_list(ob_id_src)
        else
          cells(1) = ob_id_src
        endif
        wts(1) = 1.
      elseif (ob_id_dmd > 0 .and. allocated(hru_num_cells)) then
        !active cells only (inactive and boundary cells are not pumped)
        allocate (cells(hru_num_cells(ob_id_dmd)), wts(hru_num_cells(ob_id_dmd)))
        do i = 1, hru_num_cells(ob_id_dmd)
          cell_id = hru_cells(ob_id_dmd,i)
          if (cell_id <= 0) cycle
          if (gw_state(cell_id)%stat /= 1) cycle
          ncl = ncl + 1
          cells(ncl) = cell_id
          wts(ncl) = hru_cells_fract(ob_id_dmd,i)
        enddo
        if (ncl > 0) then
          sum_wts = sum(wts(1:ncl))
          if (sum_wts > 0.) then
            wts(1:ncl) = wts(1:ncl) / sum_wts
          else
            wts(1:ncl) = 1. / real(ncl)
          endif
        endif
      endif
      if (ob_id_src > 0) then
        if (cells(1) <= 0) then
          ncl = 0
        elseif (gw_state(cells(1))%stat /= 1) then
          ncl = 0
        endif
      endif

      !no connected active cells: the whole demand is unmet (recorded for the HRU only)
      if (ncl == 0) then
        dmd_unmet = demand_vol
        call gwflow_pump_hru_sum
        return
      endif

      allocate (avail(ncl), req(ncl), removed(ncl))

      !groundwater available in each cell (storage less the pumping already taken today; ppag is negative)
      do i = 1, ncl
        cell_id = cells(i)
        avail(i) = 0.
        if (cell_id > 0) then
          if (gw_state(cell_id)%stat == 1 .and. gw_state(cell_id)%head > gw_state(cell_id)%botm) then
            avail(i) = (gw_state(cell_id)%head - gw_state(cell_id)%botm) * gw_state(cell_id)%area * &
                       gw_state(cell_id)%spyd + gw_hyd_ss(cell_id)%ppag
            avail(i) = max(0., avail(i))
          endif
        endif
      enddo

      !first split by intersected area
      do i = 1, ncl
        req(i) = demand_vol * wts(i)
        removed(i) = min(req(i), avail(i))
      enddo

      !second split: unmet demand from the cells that still have storage
      remaining = demand_vol - sum(removed)
      spare = sum(avail - removed)
      if (remaining > 0. .and. spare > 0.) then
        add = min(remaining, spare)
        do i = 1, ncl
          removed(i) = removed(i) + add * (avail(i) - removed(i)) / spare
          removed(i) = min(removed(i), avail(i))
        enddo
      endif

      extracted = sum(removed)
      dmd_unmet = max(0., demand_vol - extracted)

      !remove the water (and its heat and solutes) from each cell
      shortfall = 0.
      do i = 1, ncl
        shortfall = shortfall + max(0., req(i) - removed(i))
      enddo
      do i = 1, ncl
        cell_id = cells(i)
        if (cell_id <= 0) cycle
        call gwflow_pump_cell (cell_id, removed(i))
        !unmet demand by cell, in proportion to each cell's shortfall
        if (dmd_unmet > 0. .and. shortfall > 0.) then
          add = dmd_unmet * max(0., req(i) - removed(i)) / shortfall
          gw_hyd_ss(cell_id)%ppdf = gw_hyd_ss(cell_id)%ppdf + add
          gw_hyd_ss_mo(cell_id)%ppdf = gw_hyd_ss_mo(cell_id)%ppdf + add
          gw_hyd_ss_yr(cell_id)%ppdf = gw_hyd_ss_yr(cell_id)%ppdf + add
        endif
      enddo

      call gwflow_pump_hru_sum

      return

      contains

      !sum the demand, pumping and unmet demand for the HRU
      subroutine gwflow_pump_hru_sum
        if (ob_id_dmd > 0) then
          hru_pump(ob_id_dmd) = hru_pump(ob_id_dmd) + extracted
          hru_pump_dmd(ob_id_dmd) = hru_pump_dmd(ob_id_dmd) + demand_vol
          hru_pump_unmet(ob_id_dmd) = hru_pump_unmet(ob_id_dmd) + dmd_unmet
        endif
      end subroutine gwflow_pump_hru_sum

      !remove a volume of groundwater from one cell, with its heat and solute mass;
      !the solute mass goes to the soil (or wetland) of the HRU
      subroutine gwflow_pump_cell (icell, vol)
        integer, intent (in) :: icell       !       |gwflow cell
        real, intent (in) :: vol            !m3     |groundwater removed
        integer :: s                        !       |solute counter
        integer :: isalt                    !       |salt ion counter
        integer :: ics                      !       |constituent counter
        integer :: sol_index                !       |solute index
        integer :: wetland                  !       |wetland flag
        real :: heat_flux                   !J      |heat in pumped groundwater
        real :: gw_mass                     !kg     |solute mass left in the cell today
        real :: irr_mass(100)               !kg     |solute mass removed in the pumped water

        if (vol <= 0.) return

        !pumping volume (m3, negative = leaving the aquifer), for gwflow_simulate and output
        gw_hyd_ss(icell)%ppag = gw_hyd_ss(icell)%ppag - vol
        gw_hyd_ss_mo(icell)%ppag = gw_hyd_ss_mo(icell)%ppag - vol
        gw_hyd_ss_yr(icell)%ppag = gw_hyd_ss_yr(icell)%ppag - vol

        !heat removed with the pumped water
        if (gw_heat_flag == 1) then
          heat_flux = gwheat_state(icell)%temp * gw_rho * gw_cp * vol !J
          heat_flux = min(heat_flux, gwheat_state(icell)%stor + gw_heat_ss(icell)%ppag)
          heat_flux = max(0., heat_flux)
          gw_heat_ss(icell)%ppag = gw_heat_ss(icell)%ppag - heat_flux
          gw_heat_ss_yr(icell)%ppag = gw_heat_ss_yr(icell)%ppag - heat_flux
        endif

        !solute mass removed with the pumped water
        !(mass is removed from the aquifer via the mass balance equation in gwflow_simulate)
        if (gw_solute_flag == 1) then
          irr_mass = 0.
          do s = 1, gw_nsolute
            gw_mass = (gwsol_state(icell)%solute(s)%mass + gwsol_ss(icell)%solute(s)%ppag) / 1000. !kg
            irr_mass(s) = (gwsol_state(icell)%solute(s)%conc * vol) / 1000. !kg
            irr_mass(s) = min(irr_mass(s), gw_mass)
            irr_mass(s) = max(0., irr_mass(s))
            gwsol_ss(icell)%solute(s)%ppag = gwsol_ss(icell)%solute(s)%ppag - irr_mass(s) * 1000. !g
            gwsol_ss_sum(icell)%solute(s)%ppag = gwsol_ss_sum(icell)%solute(s)%ppag - irr_mass(s) * 1000.
            gwsol_ss_sum_mo(icell)%solute(s)%ppag = gwsol_ss_sum_mo(icell)%solute(s)%ppag - irr_mass(s) * 1000.
          enddo
          if (ob_id_dmd > 0) then
            wetland = hru(ob_id_dmd)%dbs%surf_stor
            if (wetland > 0) then !add to wetland
              wet(ob_id_dmd)%no3 = wet(ob_id_dmd)%no3 + irr_mass(1) !kg
              wet(ob_id_dmd)%solp = wet(ob_id_dmd)%solp + irr_mass(2) !kg
              sol_index = 2
              if (gwsol_salt == 1) then
                do isalt = 1, cs_db%num_salts
                  sol_index = sol_index + 1
                  wet_water(ob_id_dmd)%salt(isalt) = wet_water(ob_id_dmd)%salt(isalt) + irr_mass(sol_index) !kg
                  wetsalt_d(ob_id_dmd)%salt(isalt)%irrig = wetsalt_d(ob_id_dmd)%salt(isalt)%irrig + irr_mass(sol_index)
                enddo
              endif
              if (gwsol_cons == 1) then
                do ics = 1, cs_db%num_cs
                  sol_index = sol_index + 1
                  wet_water(ob_id_dmd)%cs(ics) = wet_water(ob_id_dmd)%cs(ics) + irr_mass(sol_index) !kg
                  wetcs_d(ob_id_dmd)%cs(ics)%irrig = wetcs_d(ob_id_dmd)%cs(ics)%irrig + irr_mass(sol_index)
                enddo
              endif
            else !add to the top soil layer (kg/ha)
              soil1(ob_id_dmd)%mn(1)%no3 = soil1(ob_id_dmd)%mn(1)%no3 + (irr_mass(1)/hru(ob_id_dmd)%area_ha)
              soil1(ob_id_dmd)%mp(1)%lab = soil1(ob_id_dmd)%mp(1)%lab + (irr_mass(2)/hru(ob_id_dmd)%area_ha)
              sol_index = 2
              if (gwsol_salt == 1) then
                do isalt = 1, cs_db%num_salts
                  sol_index = sol_index + 1
                  cs_soil(ob_id_dmd)%ly(1)%salt(isalt) = cs_soil(ob_id_dmd)%ly(1)%salt(isalt) + &
                    (irr_mass(sol_index)/hru(ob_id_dmd)%area_ha)
                  hsaltb_d(ob_id_dmd)%salt(isalt)%irgw = hsaltb_d(ob_id_dmd)%salt(isalt)%irgw + &
                    (irr_mass(sol_index)/hru(ob_id_dmd)%area_ha)
                enddo
              endif
              if (gwsol_cons == 1) then
                do ics = 1, cs_db%num_cs
                  sol_index = sol_index + 1
                  cs_soil(ob_id_dmd)%ly(1)%cs(ics) = cs_soil(ob_id_dmd)%ly(1)%cs(ics) + &
                    (irr_mass(sol_index)/hru(ob_id_dmd)%area_ha)
                  hcsb_d(ob_id_dmd)%cs(ics)%irgw = hcsb_d(ob_id_dmd)%cs(ics)%irgw + &
                    (irr_mass(sol_index)/hru(ob_id_dmd)%area_ha)
                enddo
              endif
            endif
          endif
        endif

      end subroutine gwflow_pump_cell

      end subroutine gwflow_pump_allo
