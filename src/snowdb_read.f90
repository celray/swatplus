      subroutine snowdb_read
      
      use input_file_module
      use maximum_data_module
      use hru_module
      use utils, only: open_file
      
      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: msno = 0             !           |
      integer :: isno = 0             !none       |counter
      
      msno = 0
      eof = 0
      imax = 0
      
      
      !! read snow database data from snow.sno
      if (open_file(107, in_parmdb%snow)) then
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
          
        rewind (107)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
        
        allocate (snodb(0:imax))
    
        do isno = 1, imax
          read (107,*,iostat=eof) snodb(isno)         
          if (eof < 0) exit
        end do

      exit
      enddo
      else
        allocate (snodb(0:0))
      endif
      close (107)
      
      db_mx%sno = imax
      
      return
      end subroutine snowdb_read