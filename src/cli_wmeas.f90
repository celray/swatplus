      subroutine cli_wmeas
      
      use input_file_module
      use climate_module
      use maximum_data_module
      use time_module
      use utils, only: open_file
      
      implicit none
      
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      character (len=256) :: cli_file = ""!         |station file with path prefix
      integer :: eof = 0              !           |end of file
      integer :: i = 0                !none       |counter  
      integer :: imax = 0             !none       |ending of loop
      integer :: iyr = 0              !none       |number of years 
      integer :: istep = 0            !           |
      integer :: iyr_prev = 0          !none      |previous year
      integer :: iyrs = 0             !           |

       eof = 0
       imax = 0

      !! read all measured daily wind data
      !! prepend wnd_path to station file if set
      if (trim(in_cli%wnd_cli) == "null" .or. len_trim(in_cli%wnd_cli) == 0  &
          .or. trim(in_path_wnd%wnd) == "null" .or. len_trim(in_path_wnd%wnd) == 0) then
        cli_file = in_cli%wnd_cli
      else
        cli_file = TRIM(ADJUSTL(in_path_wnd%wnd))//TRIM(in_cli%wnd_cli)
      end if
      if (open_file(107, cli_file)) then
      do
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
          do while (eof == 0)
            read (107,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do
          
      allocate (wnd(0:imax))
      allocate (wnd_n(imax))
      
      rewind (107)
      read (107,*,iostat=eof) titldum
      if (eof < 0) exit
      read (107,*,iostat=eof) header
      if (eof < 0) exit
      do i = 1, imax
        read (107,*,iostat=eof) wnd_n(i)
        if (eof < 0) exit
      end do
      
      rewind (107)
      read (107,*,iostat=eof) titldum
      if (eof < 0) exit
      read (107,*,iostat=eof) header
      if (eof < 0) exit
      
      do i = 1, imax
        read (107,*,iostat = eof) wnd(i)%filename
        if (eof < 0) exit
        
!!!!!weather path code
       if (in_path_wnd%wnd == "null") then
         open (108,file = wnd(i)%filename)
       else
        open (108,file = TRIM(ADJUSTL(in_path_wnd%wnd))//wnd(i)%filename)
       endif
!!!!!weather path code
        
        read (108,*,iostat=eof) titldum
        if (eof < 0) exit
        read (108,*,iostat=eof) header
        if (eof < 0) exit
        read (108,*,iostat=eof) wnd(i)%nbyr, wnd(i)%tstep, wnd(i)%lat, wnd(i)%long,     &
                               wnd(i)%elev
        if (eof < 0) exit
       
        ! the precip time step has to be the same as time%step
        allocate (wnd(i)%ts(366,wnd(i)%nbyr))
        
        ! read and save start jd and yr
        read (108,*,iostat=eof) iyr, istep
        if (eof < 0) exit
        
        wnd(i)%start_day = istep
        wnd(i)%start_yr = iyr
        
        backspace (108)

      if (iyr > time%yrc) then
        wnd(i)%yrs_start = iyr - time%yrc
      else
        ! read and store entire year
        wnd(i)%yrs_start = 0
      end if
      
        ! read and store entire year
       do 
         read (108,*,iostat=eof) iyr, istep
         if (eof < 0) exit
         if (iyr >= time%yrc .and. istep >= time%day_start) exit
       end do
 
       backspace (108)
       iyr_prev = iyr
       iyrs = 1
       
       do
         read (108,*,iostat=eof) iyr, istep, wnd(i)%ts(istep,iyrs)
         if (eof < 0) exit
         if (istep == 365 .or. istep == 366) then
           read (108,*,iostat=eof) iyr, istep
           if (eof < 0) exit
           backspace (108)
           if (iyr /= iyr_prev) then
             iyr_prev = iyr
             iyrs = iyrs + 1
           end if
         end if
       end do
       close (108)
      
       !save end jd and year
       wnd(i)%end_day = istep
       wnd(i)%end_yr = iyr
       
      end do
      close (107)
      exit
      end do
      else
        allocate (wnd(0:0))
        allocate (wnd_n(0))
      end if
             
       db_mx%wndfiles = imax
       
      return
      end subroutine cli_wmeas