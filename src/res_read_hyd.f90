      subroutine res_read_hyd
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use reservoir_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""  !             |title of file
      character (len=80) :: header = "" !             |header of file
      integer :: eof = 0                !             |end of file
      integer :: imax = 0               !             |determine max number for array (imax) and total number in file
      integer :: ires = 0               !none         |counter

      eof = 0
      imax = 0

      if (open_file(105, in_res%hyd_res)) then
      do
       read (105,*,iostat=eof) titldum
       if (eof < 0) exit
       read (105,*,iostat=eof) header
       if (eof < 0) exit
        do while (eof == 0)
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          imax = imax + 1
        end do
        
      db_mx%res_hyd = imax
      
      allocate (res_hyddb(0:imax))
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
      
       do ires = 1, imax
         
         !read (105,*,iostat=eof) titldum
         !backspace (105)
         read (105,*,iostat=eof) res_hyddb(ires)
         if (eof < 0) exit

        if (res_hyddb(ires)%pvol + res_hyddb(ires)%evol > 0.) then
          if(res_hyddb(ires)%pvol <= 0) res_hyddb(ires)%pvol = 0.9 * res_hyddb(ires)%evol
        else
          if (res_hyddb(ires)%pvol <= 0) res_hyddb(ires)%pvol = 60000.0
        end if
        if (res_hyddb(ires)%evol <= 0.0) res_hyddb(ires)%evol = 1.11 * res_hyddb(ires)%pvol
        if (res_hyddb(ires)%psa <= 0.0) res_hyddb(ires)%psa = 0.08 * res_hyddb(ires)%pvol
        if (res_hyddb(ires)%esa <= 0.0) res_hyddb(ires)%esa = 1.5 * res_hyddb(ires)%psa
        if (res_hyddb(ires)%evrsv <= 0.) res_hyddb(ires)%evrsv = 0.6

       end do
       close (105)
      exit
      enddo
      else
        allocate (res_hyddb(0:0))
      endif
  
      return
      end subroutine res_read_hyd