      subroutine gwflow_pump_ext !rtb gwflow

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine determines the volume of groundwater that is extracted
!!    from gwflow grid cells, for groundwater that is lost to the system
!!    (pumping volume are used in gwflow_simulate, in groundwater balance equations)

      use gwflow_module

      implicit none

      integer :: i = 0                      !       |counter
      integer :: j = 0                      !       |counter
      integer :: s = 0                      !       |solute counter
      integer :: cell_id = 0                !       |gwflow cell
      integer :: pumpex_start_date = 0      !       |starting date of pumping period
      integer :: pumpex_end_date = 0        !       |ending date of pumping period
      real :: Q = 0.                        !m3     |pumping rate
      real :: solmass(100) = 0.             !g      |solute mass in pumped groundwater
      real :: heat_flux = 0.                !J      |heat in pumped groundwater
      real :: Q_req = 0.                    !m3     |requested pumping
      real :: avail = 0.                    !m3     |groundwater left in the cell after today's pumping so far


      !only proceed if external pumping has been specified
      if (gw_pumpex_flag == 1) then

        !loop through pumps and specified pumping periods; apply rate to the associated grid cells
        do i=1,gw_npumpex

          !cell where pumping occurs
          cell_id = gw_pumpex_cell(i)

          !only proceed if pumping cell is active
          if(gw_state(cell_id)%stat == 1) then

            !loop through the pumping periods of the cell
            do j=1,gw_pumpex_nperiods(i)

              !determine if the current day of the simulation is within the pumping period; if so, apply the pumping rate to the cell
              pumpex_start_date = gw_pumpex_dates(i,1,j)
              pumpex_end_date = gw_pumpex_dates(i,2,j)
              if(gw_daycount.ge.pumpex_start_date .and. gw_daycount.le.pumpex_end_date) then

                !check to make sure there is enough groundwater to satisfy the pumping rate
                !storage less what irrigation (ppag) and earlier external pumping (ppex) took today;
                !the shortfall is recorded as unmet pumping (ppdf), as for irrigation
                Q_req = gw_pumpex_rates(i,j)
                avail = gw_state(cell_id)%stor + gw_hyd_ss(cell_id)%ppag + gw_hyd_ss(cell_id)%ppex
                avail = max(0., avail)
                Q = min(Q_req, avail)
                Q = max(0., Q)
                if(Q_req > Q) then
                  gw_hyd_ss(cell_id)%ppdf = gw_hyd_ss(cell_id)%ppdf + (Q_req - Q)
                  gw_hyd_ss_mo(cell_id)%ppdf = gw_hyd_ss_mo(cell_id)%ppdf + (Q_req - Q)
                  gw_hyd_ss_yr(cell_id)%ppdf = gw_hyd_ss_yr(cell_id)%ppdf + (Q_req - Q)
                endif
                gw_hyd_ss(cell_id)%ppex = gw_hyd_ss(cell_id)%ppex - Q !negative = leaving the aquifer
                gw_hyd_ss_yr(cell_id)%ppex = gw_hyd_ss_yr(cell_id)%ppex - Q !store for annual water
                gw_hyd_ss_mo(cell_id)%ppex = gw_hyd_ss_mo(cell_id)%ppex - Q !store for monthly water

                !heat
                if(gw_heat_flag == 1) then
                  heat_flux = gwheat_state(cell_id)%temp * gw_rho * gw_cp * Q !J
                  heat_flux = min(heat_flux, gwheat_state(cell_id)%stor + gw_heat_ss(cell_id)%ppag + &
                                  gw_heat_ss(cell_id)%ppex)
                  heat_flux = max(0., heat_flux)
                  gw_heat_ss(cell_id)%ppex = gw_heat_ss(cell_id)%ppex - heat_flux !negative = leaving the aquifer
                  gw_heat_ss_yr(cell_id)%ppex = gw_heat_ss_yr(cell_id)%ppex - heat_flux
                endif

                !if chemical transport simulated, calculate the mass of N and P removed via pumping
                if (gw_solute_flag == 1) then
                  do s=1,gw_nsolute !loop through the solutes
                    solmass(s) = Q * gwsol_state(cell_id)%solute(s)%conc !g
                    solmass(s) = min(solmass(s), gwsol_state(cell_id)%solute(s)%mass + &
                      gwsol_ss(cell_id)%solute(s)%ppag + gwsol_ss(cell_id)%solute(s)%ppex)
                    solmass(s) = max(0., solmass(s))
                    gwsol_ss(cell_id)%solute(s)%ppex = gwsol_ss(cell_id)%solute(s)%ppex - solmass(s)
                    gwsol_ss_sum(cell_id)%solute(s)%ppex = gwsol_ss_sum(cell_id)%solute(s)%ppex - solmass(s)
									gwsol_ss_sum_mo(cell_id)%solute(s)%ppex = gwsol_ss_sum_mo(cell_id)%solute(s)%ppex - solmass(s)
                  enddo
                endif

              endif

            enddo !go to next pumping period

          endif !check if cell is active

        enddo !go to next pump

      endif

      end subroutine gwflow_pump_ext
