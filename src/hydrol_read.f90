      subroutine hydrol_read
      
      use input_file_module
      use maximum_data_module
      use hydrology_data_module
      use utils, only: open_file

      implicit none

      integer :: mhydrol = 0          !           |
      integer :: ithyd = 0            !none       |counter
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file

      mhydrol = 0
      eof = 0
      imax = 0

      !! read all data from hydrol.dat
      if (open_file(107, in_hyd%hydrol_hyd)) then
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
          
          allocate (hyd_db(0:imax))
        
          rewind (107)
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
                 
          do ithyd = 1, imax
             read (107,*,iostat=eof) hyd_db(ithyd)            
             if (eof < 0) exit
          end do
          exit
        enddo
      else
        allocate (hyd_db(0:0))
      endif
      close (107)
 
      db_mx%hyd = imax
      
      return
      end subroutine hydrol_read