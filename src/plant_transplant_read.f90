      subroutine plant_transplant_read
      
      use input_file_module
      use maximum_data_module
      use plant_data_module
      use utils, only: open_file

      implicit none 

      integer :: ic = 0               !none       |counter
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      logical :: i_exist              !none       |check to determine if file exists
      
      eof = 0
      imax = 0

      if (open_file(104, in_ops%transplant_ops)) then
      do
        read (104,*,iostat=eof) titldum
        if (eof < 0) exit
        read (104,*,iostat=eof) header
        if (eof < 0) exit
          do while (eof == 0)
            read (104,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do
        allocate (transpl(0:imax))
        
        rewind (104)
        read (104,*,iostat=eof) titldum
        if (eof < 0) exit
        read (104,*,iostat=eof) header
        if (eof < 0) exit
        
        do ic = 1, imax
          read (104,*,iostat=eof) transpl(ic)
          if (eof < 0) exit
        end do
        
        exit
      enddo
      else
        allocate (transpl(0:0))
      endif

      db_mx%transplant = imax
      
      close (104)
      return
      end subroutine plant_transplant_read