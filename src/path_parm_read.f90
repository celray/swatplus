      subroutine path_parm_read
      
      use input_file_module
      use pathogen_data_module, only : path_db
      use maximum_data_module
      use utils, only: open_file
      
      implicit none

      character (len=80) :: titldum = "" !              |title of file
      character (len=80) :: header = ""  !              |header of file
      integer :: ibac = 0                !none          |counter  
      integer :: eof = 0                 !              |end of file
      integer :: imax = 0                !none          |counter 

      eof = 0
      imax = 0    
      
      !! read pathogen properties
      if (open_file(107, in_parmdb%pathcom_db)) then
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

        db_mx%path = imax
          
        allocate (path_db(0:imax))
        rewind (107)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit

        do ibac = 1, db_mx%path
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          backspace (107)
          read (107,*,iostat=eof) path_db(ibac)
          if (eof < 0) exit
        end do
        exit
      enddo
      else
         allocate (path_db(0:0))
      endif
      
      close (107)
      
      return
      end subroutine path_parm_read