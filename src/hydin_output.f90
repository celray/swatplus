      subroutine hydin_output
    
      use hydrograph_module
      use time_module
      use basin_module
#ifdef SQLITE_ENABLED
      use sqlite_output_module
#endif

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine outputs hyd variables on daily, monthly and annual time steps

      implicit none
      
      integer :: iin = 0          !none          |counter
      

      do icmd = 1, sp_ob%objs
        do iin = 1, ob(icmd)%rcv_tot
        !! daily print
         if (pco%day_print == "y" .and. pco%int_day_cur == pco%int_day) then
          if (pco%hyd%d == "y") then
             write (2560,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%obtyp_in(iin),        &
              ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_d(iin)
!            write (2560,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%obtyp_in(iin),        &
!              ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_d(iin)
              if (pco%csvout == "y") then
                write (2564,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, &
                 ob(icmd)%obtyp_in(iin), ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin),     &
                 ob(icmd)%hin_d(iin)
              end if
#ifdef SQLITE_ENABLED
              if (pco%sqliteout == "y") then
                call sqlite_insert_hyd("d", "in", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hin_d(iin))
              end if
#endif
          endif
        end if
                                                    
        ob(icmd)%hin_m(iin) = ob(icmd)%hin_m(iin) + ob(icmd)%hin_d(iin)
        ob(icmd)%hin_d(iin) = hz

        !! monthly print
        if (time%end_mo == 1) then
          if (pco%hyd%m == "y") then
            write (2561,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%obtyp_in(iin),        &
             ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_m(iin)
              if (pco%csvout == "y") then
                write (2565,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, &
                    ob(icmd)%obtyp_in(iin), ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin),  &
                    ob(icmd)%hin_m(iin)
              end if
#ifdef SQLITE_ENABLED
              if (pco%sqliteout == "y") then
                call sqlite_insert_hyd("m", "in", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hin_m(iin))
              end if
#endif
          end if
          ob(icmd)%hin_y(iin) = ob(icmd)%hin_y(iin)+ ob(icmd)%hin_m(iin)
          ob(icmd)%hin_m(iin) = hz
        endif
        
        !! yearly print
        if (time%end_yr == 1) then
          if (pco%hyd%y == "y") then
            write (2562,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name,  ob(icmd)%typ, ob(icmd)%obtyp_in(iin), &
             ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_y(iin)
            if (pco%csvout == "y") then
              write (2566,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name,  ob(icmd)%typ, ob(icmd)%num, &
              ob(icmd)%obtyp_in(iin), ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_y(iin)
            endif
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_hyd("y", "in", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hin_y(iin))
            end if
#endif
          end if
          ob(icmd)%hin_a(iin) = ob(icmd)%hin_a(iin)+ ob(icmd)%hin_y(iin)
          ob(icmd)%hin_y(iin) = hz
        endif
        
        !! average annual print
        if (time%end_sim == 1 .and. pco%hyd%a == "y") then
          ob(icmd)%hin_a(iin) = ob(icmd)%hin_a(iin) / time%yrs_prt
          write (2563,*) time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name,  ob(icmd)%typ,  ob(icmd)%obtyp_in(iin), &
             ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_a(iin)
            if (pco%csvout == "y") then
              write (2567,'(*(G0.6,:","))') time%day, time%mo, time%day_mo, time%yrc, ob(icmd)%name, ob(icmd)%typ, &
                ob(icmd)%obtyp_in(iin), ob(icmd)%obtypno_in(iin), ob(icmd)%htyp_in(iin), ob(icmd)%frac_in(iin), ob(icmd)%hin_a(iin)
            end if
#ifdef SQLITE_ENABLED
            if (pco%sqliteout == "y") then
              call sqlite_insert_hyd("a", "in", time%day, time%mo, time%day_mo, time%yrc, icmd, ob(icmd)%gis_id, ob(icmd)%name, ob(icmd)%typ, ob(icmd)%hin_a(iin))
            end if
#endif
        end if

        end do
      end do
      end subroutine hydin_output