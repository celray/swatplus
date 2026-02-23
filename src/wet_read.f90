      subroutine wet_read
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use reservoir_data_module
      use reservoir_module
      use hydrograph_module
      use constituent_mass_module
      use pesticide_data_module
      use res_salt_module
      use res_cs_module
      use utils, only: open_file

      implicit none

      
      
      
      external :: search
      character (len=80) :: titldum = "" !           |title of file
      character (len=80) :: header = ""  !           |header of file
      integer :: eof = 0                 !           |end of file
      integer :: imax = 0                !none       |determine max number for array (imax) and total number in file
      integer :: i = 0                   !none       |counter
      integer :: ires = 0                !none       |counter 
      integer :: k = 0                   !none       |counter 
      
      eof = 0
      imax = 0
            
      !read wetland.wet
      imax = 0
      if (open_file(105, in_res%wet)) then
        do
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          read (105,*,iostat=eof) header
          if (eof < 0) exit
          do while (eof == 0)
            read (105,*,iostat=eof) i
            if (eof < 0) exit
            imax = imax + 1
          end do
        
          db_mx%wet_dat = imax
       
          allocate (wet_dat_c(imax))
          allocate (wet_dat(imax))
      
          rewind (105)
          read (105,*,iostat = eof) titldum
          if (eof < 0) exit
          read (105,*,iostat=eof) header
          if (eof < 0) exit
      
          do ires = 1, db_mx%wet_dat
            read (105,*,iostat=eof) i
            if (eof < 0) exit
            backspace (105)
            read (105,*,iostat=eof) k, wet_dat_c(ires)
            if (eof < 0) exit
          end do
      
          db_mx%wet_dat = imax
            
         close (105)
         exit
       end do
      else
        allocate (wet_dat_c(0:0))
        allocate (wet_dat(0:0))
      end if
      
      return
      end subroutine wet_read