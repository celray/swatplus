      subroutine septic_parm_read
      
      use input_file_module
      use maximum_data_module
      use septic_data_module
      use utils, only: open_file
      
      implicit none
         
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: is = 0               !none       |counter
      
      eof = 0
      imax = 0

      if (open_file(171, in_parmdb%septic_sep)) then
      do
        read (171,*,iostat=eof) titldum
        if (eof < 0) exit
        read (171,*,iostat=eof) header
        if (eof < 0) exit
           do while (eof == 0) 
             read (171,*,iostat=eof) titldum
             if (eof < 0) exit
             imax = imax + 1
           end do
           
        db_mx%sep = imax
           
        allocate (sepdb(0:imax))
        rewind (171)
        read (171,*,iostat=eof) titldum
        if (eof < 0) exit
        read (171,*,iostat=eof) header
        if (eof < 0) exit
    
        do is = 1, db_mx%sep
          read (171,*,iostat=eof) titldum
          if (eof < 0) exit
          backspace (171)
          read (171,*,iostat=eof) sepdb(is)
          if (eof < 0) exit
        end do

        close (171)
        exit
      enddo
      else
          allocate (sepdb(0:0))
      endif
      
      return
      end subroutine septic_parm_read