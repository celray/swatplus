      subroutine aquifer_output(iaq)
      
      use time_module
      use basin_module
      use aquifer_module
      use hydrograph_module, only : ob, sp_ob1
#ifdef SQLITE_ENABLED
      use sqlite_output_module
#endif
      
      implicit none
            
      integer, intent (in) :: iaq        !             |
      real :: const = 0.                 !             |constant used for rate, days, etc
      integer :: iob = 0                 !             |
                          
      iob = sp_ob1%aqu + iaq - 1
       
        !! sum monthly variables        
        aqu_m(iaq) = aqu_m(iaq) + aqu_d(iaq)
        
        !! daily print - AQUIFER
         if (pco%day_print == "y" .and. pco%int_day_cur == pco%int_day) then
          if (pco%aqu%d == "y") then
            if (pco%textout == "y") then
              write (2520,100) time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_d(iaq)
              if (pco%csvout == "y") then
                write (2524,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_d(iaq)
              end if
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_aquifer("d", time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_d(iaq))
            end if
#endif
          end if
        end if

        !! monthly print - AQUIFER
        if (time%end_mo == 1) then
          const = float (ndays(time%mo + 1) - ndays(time%mo))
          aqu_m(iaq)%stor = aqu_m(iaq)%stor / const 
          aqu_m(iaq)%dep_wt = aqu_m(iaq)%dep_wt / const
          aqu_m(iaq)%no3_st = aqu_m(iaq)%no3_st / const
          aqu_y(iaq) = aqu_y(iaq) + aqu_m(iaq)
          if (pco%aqu%m == "y") then
            if (pco%textout == "y") then
              write (2521,100)  time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_m(iaq)
              if (pco%csvout == "y") then
                write (2525,'(*(G0.6,:","))')  time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_m(iaq)
              endif
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_aquifer("m", time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_m(iaq))
            end if
#endif
          end if
          aqu_m(iaq) = aquz
        end if

        !! yearly print - AQUIFER
        if (time%end_yr == 1) then
          aqu_y(iaq)%stor = aqu_y(iaq)%stor / 12.
          aqu_y(iaq)%dep_wt = aqu_y(iaq)%dep_wt / 12.
          aqu_y(iaq)%no3_st = aqu_y(iaq)%no3_st / 12.
          aqu_a(iaq) = aqu_a(iaq) + aqu_y(iaq)
          if (pco%aqu%y == "y") then
            if (pco%textout == "y") then
              write (2522,102) time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_y(iaq)
              if (pco%csvout == "y") then
                write (2526,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_y(iaq)
              end if
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_aquifer("y", time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_y(iaq))
            end if
#endif
          end if
          !! zero yearly variables        
          aqu_y(iaq) = aquz
        end if
        
      !! average annual print - AQUIFER
      if (time%end_sim == 1 .and. pco%aqu%a == "y") then
        aqu_a(iaq) = aqu_a(iaq) / time%yrs_prt
        if (pco%textout == "y") then
          write (2523,102) time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_a(iaq)
          if (pco%csvout == "y") then 
            write (2527,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_a(iaq)
          end if 
        end if
#ifdef SQLITE_ENABLED
        if (pco%sqliteout == "y") then
          call sqlite_insert_aquifer("a", time%day, time%mo, time%day_mo, time%yrc, iaq, ob(iob)%gis_id, ob(iob)%name, aqu_a(iaq))
        end if
#endif
      end if
      
      return
      
100   format (4i6,2i8,2x,a,20f15.3)
102   format (4i6,2i8,2x,a,20f15.3)
       
      end subroutine aquifer_output