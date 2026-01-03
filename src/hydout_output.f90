      subroutine hydout_output (iout)
    
      use time_module
      use basin_module
      use hydrograph_module
#ifdef SQLITE_ENABLED
      use sqlite_output_module
#endif
      
      implicit none
      
      integer :: iout           !           |

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine outputs hyd variables on daily, monthly and annual time steps
      
      !!  0 = average annual (always print)
      !!  1 = yearly
      !!  2 = monthly
      !!  3 = daily  

!!!!! daily print
         if (pco%day_print == "y" .and. pco%int_day_cur == pco%int_day) then
          if (pco%hyd%d == "y") then
            write (2580,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ,          &
              ob(icmd)%obtyp_out(iout),                      &
             ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),           &
             ob(icmd)%frac_out(iout), ht1
            if (pco%csvout == "y") then
              write (2584,'(*(G0.3,:","))')time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ,   &
               ob(icmd)%obtyp_out(iout),                              &
               ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),                   &
               ob(icmd)%frac_out(iout), ht1  
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_hyd("d", "out", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ht1)
            end if
#endif
          end if
        endif
        ob(icmd)%hout_m(iout) = ob(icmd)%hout_m(iout) + ht1

!!!!! monthly print
        if (time%end_mo == 1) then
          if (pco%hyd%m == "y") then
            write (2581,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ,      & 
           ob(icmd)%obtyp_out(iout),                      &
           ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),           &
           ob(icmd)%frac_out(iout), ob(icmd)%hout_m(iout)
            if (pco%csvout == "y") then
              write (2585,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, & 
             ob(icmd)%obtyp_out(iout),                              &
             ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),                   &
             ob(icmd)%frac_out(iout), ob(icmd)%hout_m(iout)
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_hyd("m", "out", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hout_m(iout))
            end if
#endif
          end if
          ob(icmd)%hout_y(iout) = ob(icmd)%hout_y(iout) +                 &
                                             ob(icmd)%hout_m(iout)
          ob(icmd)%hout_m(iout) = hz
        end if
        
!!!!! yearly print
        if (time%end_yr == 1) then
          if (pco%hyd%y == "y") then
            write (2582,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ,        &
            ob(icmd)%obtyp_out(iout),                      &
           ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),           &
           ob(icmd)%frac_out(iout), ob(icmd)%hout_y(iout)
             if (pco%csvout == "y") then
               write (2586,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ,  &
               ob(icmd)%obtyp_out(iout),                              &
               ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),                   &
               ob(icmd)%frac_out(iout), ob(icmd)%hout_y(iout)
             end if
#ifdef SQLITE_ENABLED
             if (pco%sqliteout == "y") then
               call sqlite_insert_hyd("y", "out", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hout_y(iout))
             end if
#endif
          end if
          ob(icmd)%hout_a(iout) = ob(icmd)%hout_a(iout) + ob(icmd)%hout_y(iout)
          ob(icmd)%hout_y(iout) = hz
        end if
        
!!!!! average annual print
        if (time%end_sim == 1 .and. pco%hyd%a == "y") then
          ob(icmd)%hout_a(iout) = ob(icmd)%hout_a(iout) / time%yrs_prt
          write (2583,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name,       &
           ob(icmd)%typ, ob(icmd)%obtyp_out(iout),        &
           ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),           &
           ob(icmd)%frac_out(iout), ob(icmd)%hout_a(iout)
            if (pco%csvout == "y") then
              write (2587,'(*(G0.3,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name,    &
              ob(icmd)%typ, ob(icmd)%obtyp_out(iout),                   &
              ob(icmd)%obtypno_out(iout), ob(icmd)%htyp_out(iout),                      &
              ob(icmd)%frac_out(iout), ob(icmd)%hout_a(iout)
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_hyd("a", "out", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hout_a(iout))
            end if
#endif
        end if
        
      return

      end subroutine hydout_output