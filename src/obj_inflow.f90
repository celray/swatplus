      subroutine obj_inflow (iob)
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine moves the current day of the object hydrographs and sums
!!    the hydrographs the object (icmd) receives
!!    (moved out of command so batches of hrus can call it in parallel)

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use time_module
      use hydrograph_module
      use constituent_mass_module
      use recall_module
      implicit none

      integer, intent (inout) :: iob  !              |last incoming object (left unchanged if there is none)
      real, dimension(time%step) :: hyd_flo     !flow hydrograph
      integer :: in                   !              | 
      integer :: iday                 !              |
      integer :: irec                 !              |
      integer :: ihyd                 !              |
      real :: conv                    !              |
      real :: frac_in                 !              |
      real :: sumflo

        if (ob(icmd)%typ == "hru" .or. ob(icmd)%typ == "ru") then
          !! hru and ru can have hyrdographs that lag into next day
          ob(icmd)%day_cur = ob(icmd)%day_cur + 1
          if (ob(icmd)%day_cur > ob(icmd)%day_max) ob(icmd)%day_cur = 1
        else
        !! assume only one day is saved for all other objects
          ob(icmd)%day_cur = 1
          !!update current day of hydrograph for the object
          ob(icmd)%day_cur = ob(icmd)%day_cur + 1
          if (ob(icmd)%day_cur > ob(icmd)%day_max) ob(icmd)%day_cur = 1
        end if
        
        !sum all receiving hydrographs
        !if (ob(icmd)%rcv_tot > 0) then
          ob(icmd)%hin = hz
          ob(icmd)%hin_sur = hz
          ob(icmd)%hin_lat = hz
          ob(icmd)%hin_til = hz
          ht1 = hz
          if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
            obcs(icmd)%hin = hin_csz
            obcs(icmd)%hin_sur = hin_csz
            obcs(icmd)%hin_lat = hin_csz
            obcs(icmd)%hin_til = hin_csz
          endif
          hcs1 = hin_csz
          hcs2 = hin_csz
          hcs3 = hin_csz
          ob(icmd)%tsin = 0.
          ob(icmd)%peakrate = 0.
          hyd_flo = 0.
          
          if (ob(icmd)%rcv_tot > 0) then
          do in = 1, ob(icmd)%rcv_tot
            iob = ob(icmd)%obj_in(in)
            ihyd = ob(icmd)%ihtyp_in(in)
            frac_in = ob(icmd)%frac_in(in)
            ob(icmd)%peakrate = ob(iob)%peakrate
            
            ! if object is not an hru, need ht1, don't need %hin_sur and %hin_lat
            ! don't have to check if it's in an ru - only hru's can be routed over
            if (ob(icmd)%typ == "hru" .or. ob(icmd)%typ == "ru" .or. ob(icmd)%typ == "hru_lte") then
                
              !! if incoming object is not an hru or ru, send total hyd to surface runoff
              if (ob(icmd)%obtyp_in(in) == "hru" .or. ob(icmd)%obtyp_in(in) == "ru" .or.          &
                                                       ob(icmd)%obtyp_in(in) == "hru_lte") then
                ! receiving hru, needs %hin_sur and %hin_lat and %hin_til to route separately in hru_control
                if (ob(icmd)%htyp_in(in) == "tot") then
                  ! if total hyd coming in from hru or ru -> add both surface and lateral flows
                  ! add to surface runon
                  ob(icmd)%hin_sur = ob(icmd)%hin_sur + frac_in * ob(iob)%hd(3)
                  if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                    obcs(icmd)%hin_sur(1) = obcs(icmd)%hin_sur(1) + frac_in * obcs(iob)%hd(3)
                  end if
                  ! add to tile flow
                  ob(icmd)%hin_til = ob(icmd)%hin_til + frac_in * ob(iob)%hd(5)
                  if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                    obcs(icmd)%hin_til(1) = obcs(icmd)%hin_til(1) + frac_in * obcs(iob)%hd(5)
                  end if
                  ! add to lateral soil runon
                  ob(icmd)%hin_lat = ob(icmd)%hin_lat + frac_in * ob(iob)%hd(4)
                  if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                    obcs(icmd)%hin_lat(1) = obcs(icmd)%hin_lat(1) + frac_in * obcs(iob)%hd(4)
                  end if
                else
                  ! if hyd in is not a total hyd from an hru or ru -> add the specified hyd typ 
                  select case (ob(icmd)%htyp_in(in))
                  case ("tot")   ! total flow
                    ob(icmd)%hin_sur = ob(icmd)%hin_sur + frac_in * ob(iob)%hd(ihyd)
                    !add constituents
                    if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                      obcs(icmd)%hin_til(1) = obcs(icmd)%hin_til(1) + frac_in * obcs(iob)%hd(ihyd)
                    end if
                  case ("sur")   ! surface runoff
                    ob(icmd)%hin_sur = ob(icmd)%hin_sur + frac_in * ob(iob)%hd(ihyd)
                    !add constituents
                    if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                      obcs(icmd)%hin_sur(1) = obcs(icmd)%hin_sur(1) + frac_in * obcs(iob)%hd(ihyd)
                    end if
                  case ("lat")   ! lateral soil flow
                    ob(icmd)%hin_lat = ob(icmd)%hin_lat + frac_in * ob(iob)%hd(ihyd)
                    !add constituents
                    if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                      obcs(icmd)%hin_lat(1) = obcs(icmd)%hin_lat(1) + frac_in * obcs(iob)%hd(ihyd)
                    end if
                  case ("til")   ! tile flow
                    ob(icmd)%hin_til = ob(icmd)%hin_til + frac_in * ob(iob)%hd(ihyd)
                    !add constituents
                    if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                      obcs(icmd)%hin_til(1) = obcs(icmd)%hin_til(1) + frac_in * obcs(iob)%hd(ihyd)
                    end if
                  case ("aqu")   ! aquifer inflow
                    ob(icmd)%hin_aqu = ob(icmd)%hin_aqu + frac_in * ob(iob)%hd(ihyd)
                    !add constituents
                    if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                      obcs(icmd)%hin_aqu(1) = obcs(icmd)%hin_aqu(1) + frac_in * obcs(iob)%hd(ihyd)
                    end if
                  end select
                end if  
              else
                ! add total inflow to surface runon if channel or recall
                ob(icmd)%hin_sur = ob(icmd)%hin_sur + frac_in * ob(iob)%hd(1)
                if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                  obcs(icmd)%hin_sur(1) = obcs(icmd)%hin_sur(1) + frac_in * obcs(iob)%hd(1)
                end if
              end if
              
            else
              ! all objects other than hru's
              ! fraction of organics
              ht1 = frac_in * ob(iob)%hd(ihyd)
              ob(icmd)%hin = ob(icmd)%hin + ht1

              !rtb hydrograph separation
              hdsep1%flo_surq = frac_in * (ob(iob)%hdsep%flo_surq)
              hdsep1%flo_latq = frac_in * (ob(iob)%hdsep%flo_latq)
              hdsep1%flo_gwsw = frac_in * (ob(iob)%hdsep%flo_gwsw)
              hdsep1%flo_swgw = frac_in * (ob(iob)%hdsep%flo_swgw)
              hdsep1%flo_satex = frac_in * (ob(iob)%hdsep%flo_satex)
              hdsep1%flo_satexsw = frac_in * (ob(iob)%hdsep%flo_satexsw)
              hdsep1%flo_tile = frac_in * (ob(iob)%hdsep%flo_tile)
              ob(icmd)%hdsep_in%flo_surq = ob(icmd)%hdsep_in%flo_surq + hdsep1%flo_surq
              ob(icmd)%hdsep_in%flo_latq = ob(icmd)%hdsep_in%flo_latq + hdsep1%flo_latq
              ob(icmd)%hdsep_in%flo_gwsw = ob(icmd)%hdsep_in%flo_gwsw + hdsep1%flo_gwsw
              ob(icmd)%hdsep_in%flo_swgw = ob(icmd)%hdsep_in%flo_swgw + hdsep1%flo_swgw
              ob(icmd)%hdsep_in%flo_satex = ob(icmd)%hdsep_in%flo_satex + hdsep1%flo_satex
              ob(icmd)%hdsep_in%flo_satexsw = ob(icmd)%hdsep_in%flo_satexsw + hdsep1%flo_satexsw
              ob(icmd)%hdsep_in%flo_tile = ob(icmd)%hdsep_in%flo_tile + hdsep1%flo_tile
              
              ! fraction of constituents
              if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                hcs1 = frac_in * obcs(iob)%hd(ihyd)
                obcs(icmd)%hin(1) = obcs(icmd)%hin(1) + hcs1
              end if
              ob(icmd)%hin_d(in) = ht1        !for hydrograph output
              if (cs_db%num_tot > 0 .and. obcs_alloc(icmd).eq.1) then
                obcs(icmd)%hcsin_d(in) = hcs1   !for constituent hydrograph output
              endif
            end if
            
            !sum subdaily inflow hydrographs
            !if (time%step > 0) then
              iday = ob(iob)%day_cur
              if (ob(iob)%typ == "hru" .or. ob(iob)%typ == "ru") then
                select case (ob(icmd)%htyp_in(in))
                case ("tot")   ! total flow
                  hyd_flo = ob(iob)%hyd_flo(iday,:) + (ob(iob)%hd(4)%flo + ob(iob)%hd(5)%flo) / time%step
                case ("sur")   ! surface runoff
                  hyd_flo(:) = ob(iob)%hyd_flo(iday,:)
                case ("rhg")   ! recharge
                  hyd_flo(:) = ob(iob)%hd(2)%flo / time%step
                case ("lat")   ! lateral soil flow
                  hyd_flo(:) = ob(iob)%hd(4)%flo / time%step
                case ("til")   ! tile flow
                  hyd_flo(:) = ob(iob)%hd(5)%flo / time%step
                end select
              end if
              select case (ob(iob)%typ)
              case ("aqu")      ! aquifer inflow
                hyd_flo(:) = ob(iob)%hd(ihyd)%flo / time%step
              case ("chandeg")  ! channel inflow
                hyd_flo(:) = ob(iob)%hyd_flo(1,:)
                sumflo = sum (hyd_flo(:))
                sumflo = 1. * sumflo
              case ("res")      ! reservoir inflow
                hyd_flo(:) = ob(iob)%hd(ihyd)%flo / time%step
              case ("outlet")      ! outlet inflow
                hyd_flo(:) = ob(iob)%hd(ihyd)%flo / time%step
              case ("recall")   ! point source inflow
                irec = ob(iob)%num
                if (recall(irec)%typ == 0) then    !subdaily
                  hyd_flo(:) = ob(iob)%hyd_flo(ob(iob)%day_cur,:)
                else                                ! monthly, yearly, and ave annual
                  hyd_flo(:) = ob(iob)%hd(1)%flo / time%step
                end if
              end select
                
              !! multiply inflow hyd by the fraction of incoming
              hyd_flo = frac_in * hyd_flo
              !! add flow hydrographs for each incoming object
              ob(icmd)%tsin = ob(icmd)%tsin + hyd_flo
              
            !end if

          end do    ! in = 1, ob(icmd)%rcv_tot

          !convert to per area basis
          if (ob(icmd)%typ == "hru" .or. ob(icmd)%typ == "ru") then  !only convert hru and subbasin hyds for routing
            conv = ob(icmd)%area_ha
            ob(icmd)%hin_sur = ob(icmd)%hin_sur / conv
            ob(icmd)%hin_sur%flo = ob(icmd)%hin_sur%flo / 10.      ! m3/10*ha = mm
            ob(icmd)%hin_lat = ob(icmd)%hin_lat / conv
            ob(icmd)%hin_lat%flo = ob(icmd)%hin_lat%flo / 10.      ! m3/10*ha = mm
            ob(icmd)%hin_til = ob(icmd)%hin_til / conv
            ob(icmd)%hin_til%flo = ob(icmd)%hin_til%flo / 10.      ! m3/10*ha = mm
          end if
        end if

      return
      end subroutine obj_inflow
