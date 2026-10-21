      module gwflow_defer_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    keeps what gwflow_soil and gwflow_wetland add to the cell sums, one log per hru,
!!    and adds it to the sums in command order (gwdef_apply, called from hru_sums).
!!    the hrus of a batch then do not touch the shared cell sums, and each sum gets the
!!    same terms in the same order as in the serial loop, so the results are identical.
!!
!!    what is logged, per hru and connected cell, in the order the serial code adds it:
!!      gwflow_soil     gw_hyd_ss%soil, gw_hyd_ss_yr%soil, gw_hyd_ss_mo%soil (water)
!!                      gw_heat_ss%soil, gw_heat_ss_yr%soil (heat, if gw_heat_flag = 1)
!!      gwflow_wetland  gw_hyd_ss%wetl, gw_hyd_ss_yr%wetl, gw_hyd_ss_mo%wetl (water)
!!                      gw_heat_ss%wetl, gw_heat_ss_yr%wetl (heat, if gw_heat_flag = 1)
!!                      gwsol_ss%wetl (set, the last hru in command order wins),
!!                      gwsol_ss_sum%wetl, gwsol_ss_sum_mo%wetl (solutes, if gw_solute_flag = 1)
!!    every logged term is the value the serial statement adds. none of the terms depends
!!    on a cell sum, or on cell state that an hru changes (gw_state, gwheat_state and
!!    gwsol_state change only in gwflow_simulate), so no gate is needed.
!!
!!    storage: slots of hru j are gwdef_off(j)+1 ... gwdef_off(j+1), one per connected cell
!!    (hru_num_cells), allocated once in gwdef_init. an hru writes only the slots of that hru.
!!    without gwdef_on (gwdef_init not called, or gwflow not active) the terms are added
!!    at once, as in the original code.

      implicit none

      logical :: gwdef_on = .false.                                   !none   |terms are logged
      integer, dimension (:), allocatable :: gwdef_off                !none   |slots before the first slot of the hru
      integer, dimension (:), allocatable :: gwdef_ns                 !none   |soil entries logged by the hru
      integer, dimension (:), allocatable :: gwdef_nw                 !none   |wetland entries logged by the hru

      !! gwflow_soil
      integer, dimension (:), allocatable :: gws_cell                 !none   |cell
      real, dimension (:), allocatable :: gws_q                       !m3     |water term (hru_Q * (-1))
      real, dimension (:), allocatable :: gws_h                       !J      |heat term (heat_flux * (-1))
      logical, dimension (:), allocatable :: gws_hflag                !none   |heat term is added

      !! gwflow_wetland
      integer, dimension (:), allocatable :: gww_cell                 !none   |cell
      real, dimension (:), allocatable :: gww_q                       !m3     |water term (gw_inflow * (-1))
      real, dimension (:), allocatable :: gww_h                       !J      |heat term (heat_flux * (-1))
      logical, dimension (:), allocatable :: gww_hflag                !none   |heat term is added
      logical, dimension (:), allocatable :: gww_sflag                !none   |solute terms are used
      real, dimension (:,:), allocatable :: gww_sol                   !g      |solute terms (solute, slot)

      contains

!!    allocate the logs; call once after gwflow_read (from batch_init)
      subroutine gwdef_init

      use basin_module, only : bsn_cc
      use hydrograph_module, only : sp_ob
      use gwflow_module, only : hru_num_cells, gw_soil_flag, gw_wet_flag, gw_solute_flag, gw_nsolute

      integer :: j                    !none   |hru
      integer :: ntot                 !none   |connections of all hrus
      integer :: nsol                 !none   |solutes

      gwdef_on = .false.
      if (bsn_cc%gwflow /= 1) return
      if (.not. allocated(hru_num_cells)) return
      if (gw_soil_flag /= 1 .and. gw_wet_flag /= 1) return

      allocate (gwdef_off(sp_ob%hru + 1), source = 0)
      allocate (gwdef_ns(sp_ob%hru), source = 0)
      allocate (gwdef_nw(sp_ob%hru), source = 0)
      do j = 1, sp_ob%hru
        gwdef_off(j+1) = gwdef_off(j) + hru_num_cells(j)
      end do
      ntot = gwdef_off(sp_ob%hru + 1)

      if (gw_soil_flag == 1) then
        allocate (gws_cell(ntot), source = 0)
        allocate (gws_q(ntot), source = 0.)
        allocate (gws_h(ntot), source = 0.)
        allocate (gws_hflag(ntot), source = .false.)
      end if
      if (gw_wet_flag == 1) then
        allocate (gww_cell(ntot), source = 0)
        allocate (gww_q(ntot), source = 0.)
        allocate (gww_h(ntot), source = 0.)
        allocate (gww_hflag(ntot), source = .false.)
        allocate (gww_sflag(ntot), source = .false.)
        nsol = 0
        if (gw_solute_flag == 1) nsol = gw_nsolute
        allocate (gww_sol(max(nsol,1), ntot), source = 0.)
      end if
      gwdef_on = .true.
      return
      end subroutine gwdef_init

!!    gwflow_soil: water and heat terms of one connected cell
      subroutine gwdef_soil (j, cell, q, hflag, h)

      integer, intent (in) :: j       !none   |hru
      integer, intent (in) :: cell    !none   |cell
      real, intent (in) :: q          !m3     |water term
      logical, intent (in) :: hflag   !none   |heat term is added
      real, intent (in) :: h          !J      |heat term
      integer :: n                    !none   |entries of the hru
      integer :: i                    !none   |slot

      if (.not. gwdef_on) then
        call gwdef_add_soil (cell, q, hflag, h)
        return
      end if
      n = gwdef_ns(j) + 1
      if (gwdef_off(j) + n > gwdef_off(j+1)) error stop "gwflow_defer: more soil entries than cells for an hru"
      gwdef_ns(j) = n
      i = gwdef_off(j) + n
      gws_cell(i) = cell
      gws_q(i) = q
      gws_h(i) = h
      gws_hflag(i) = hflag
      return
      end subroutine gwdef_soil

!!    gwflow_wetland: water, heat and solute terms of one connected cell
      subroutine gwdef_wetl (j, cell, q, hflag, h, sflag, nsol, sol)

      integer, intent (in) :: j       !none   |hru
      integer, intent (in) :: cell    !none   |cell
      real, intent (in) :: q          !m3     |water term
      logical, intent (in) :: hflag   !none   |heat term is added
      real, intent (in) :: h          !J      |heat term
      logical, intent (in) :: sflag   !none   |solute terms are used
      integer, intent (in) :: nsol    !none   |solutes
      real, dimension (:), intent (in) :: sol   !g   |solute terms
      integer :: n                    !none   |entries of the hru
      integer :: i                    !none   |slot

      if (.not. gwdef_on) then
        call gwdef_add_wetl (cell, q, hflag, h, sflag, nsol, sol)
        return
      end if
      n = gwdef_nw(j) + 1
      if (gwdef_off(j) + n > gwdef_off(j+1)) error stop "gwflow_defer: more wetland entries than cells for an hru"
      gwdef_nw(j) = n
      i = gwdef_off(j) + n
      gww_cell(i) = cell
      gww_q(i) = q
      gww_h(i) = h
      gww_hflag(i) = hflag
      gww_sflag(i) = sflag
      if (sflag) gww_sol(1:nsol, i) = sol(1:nsol)
      return
      end subroutine gwdef_wetl

!!    add the logged terms of hru j to the cell sums; called from hru_sums in command order
      subroutine gwdef_apply (j)

      use gwflow_module, only : gw_nsolute

      integer, intent (in) :: j       !none   |hru
      integer :: i                    !none   |slot

      if (.not. gwdef_on) return
      do i = gwdef_off(j) + 1, gwdef_off(j) + gwdef_ns(j)
        call gwdef_add_soil (gws_cell(i), gws_q(i), gws_hflag(i), gws_h(i))
      end do
      gwdef_ns(j) = 0
      do i = gwdef_off(j) + 1, gwdef_off(j) + gwdef_nw(j)
        call gwdef_add_wetl (gww_cell(i), gww_q(i), gww_hflag(i), gww_h(i), gww_sflag(i), gw_nsolute, gww_sol(:,i))
      end do
      gwdef_nw(j) = 0
      return
      end subroutine gwdef_apply

!!    the additions of gwflow_soil (the statements of the original code)
      subroutine gwdef_add_soil (cell, q, hflag, h)

      use gwflow_module, only : gw_hyd_ss, gw_hyd_ss_yr, gw_hyd_ss_mo, gw_heat_ss, gw_heat_ss_yr

      integer, intent (in) :: cell
      real, intent (in) :: q
      logical, intent (in) :: hflag
      real, intent (in) :: h

      gw_hyd_ss(cell)%soil = gw_hyd_ss(cell)%soil + q
      gw_hyd_ss_yr(cell)%soil = gw_hyd_ss_yr(cell)%soil + q
      gw_hyd_ss_mo(cell)%soil = gw_hyd_ss_mo(cell)%soil + q
      if (hflag) then
        gw_heat_ss(cell)%soil = gw_heat_ss(cell)%soil + h
        gw_heat_ss_yr(cell)%soil = gw_heat_ss_yr(cell)%soil + h
      end if
      return
      end subroutine gwdef_add_soil

!!    the additions of gwflow_wetland (the statements of the original code)
      subroutine gwdef_add_wetl (cell, q, hflag, h, sflag, nsol, sol)

      use gwflow_module, only : gw_hyd_ss, gw_hyd_ss_yr, gw_hyd_ss_mo, gw_heat_ss, gw_heat_ss_yr,   &
                                gwsol_ss, gwsol_ss_sum, gwsol_ss_sum_mo

      integer, intent (in) :: cell
      real, intent (in) :: q
      logical, intent (in) :: hflag
      real, intent (in) :: h
      logical, intent (in) :: sflag
      integer, intent (in) :: nsol
      real, dimension (:), intent (in) :: sol
      integer :: s                    !none   |solute counter

      gw_hyd_ss(cell)%wetl = gw_hyd_ss(cell)%wetl + q
      gw_hyd_ss_yr(cell)%wetl = gw_hyd_ss_yr(cell)%wetl + q
      gw_hyd_ss_mo(cell)%wetl = gw_hyd_ss_mo(cell)%wetl + q
      if (hflag) then
        gw_heat_ss(cell)%wetl = gw_heat_ss(cell)%wetl + h
        gw_heat_ss_yr(cell)%wetl = gw_heat_ss_yr(cell)%wetl + h
      end if
      if (sflag) then
        do s = 1, nsol
          gwsol_ss(cell)%solute(s)%wetl = sol(s)
          gwsol_ss_sum(cell)%solute(s)%wetl = gwsol_ss_sum(cell)%solute(s)%wetl + sol(s)
          gwsol_ss_sum_mo(cell)%solute(s)%wetl = gwsol_ss_sum_mo(cell)%solute(s)%wetl + sol(s)
        end do
      end if
      return
      end subroutine gwdef_add_wetl

      end module gwflow_defer_module
