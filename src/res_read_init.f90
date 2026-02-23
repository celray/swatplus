      subroutine res_read_init
      
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

      !read init
      if (open_file(105, in_res%init_res)) then
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
        
      db_mx%res_init = imax
      
      allocate (res_init(0:imax))
      allocate (wet_init(0:imax))
      allocate (res_init_dat_c(0:imax))
      
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
           
       do ires = 1, db_mx%res_init
         read (105,*,iostat=eof) res_init_dat_c(ires)
         if (eof < 0) exit
       end do
       close (105)
      exit
      enddo
      else
        allocate (res_init(0:0))
        allocate (wet_init(0:0))
      endif
      
      return
      end subroutine res_read_init