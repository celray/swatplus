      subroutine solt_db_read
      
      use input_file_module
      use maximum_data_module
      use soil_data_module
      use utils, only: open_file
      
      implicit none
      
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: msolt_db = 0         !           |
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: isolt = 0            !           |
      
      msolt_db = 0
      eof = 0
      imax = 0
      
      !! read all soil test operations data from soiltest.dat
      if (open_file(107, in_sol%nut_sol)) then
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
                
          allocate (solt_db(0:imax))
         
          rewind (107)
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
                
          do isolt = 1, imax
            read (107,*,iostat=eof) solt_db(isolt)
            if (eof < 0) exit
            if (solt_db(isolt)%exp_co > 0.005) solt_db(isolt)%exp_co = 0.001
          end do 
          exit
        enddo
      else
        allocate (solt_db(0:0))
      endif
      
      db_mx%soiltest = imax
      
      close(107)      
      return  
      end subroutine solt_db_read