      subroutine sep_read

      use input_file_module
      use maximum_data_module
      use septic_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: isep = 0             !none       |counter

      eof = 0
      imax = 0

      if (open_file(172, in_str%septic_str)) then
        do
          read (172,*,iostat=eof) titldum
          if (eof < 0) exit
          read (172,*,iostat=eof) header
          if (eof < 0) exit
          do while (eof == 0)
            read (172,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do

          allocate (sep(0:imax))
          rewind (172)
          read (172,*,iostat=eof) titldum
          if (eof < 0) exit
          read (172,*,iostat=eof) header
          if (eof < 0) exit

          do isep = 1, imax
            read(172,*,iostat=eof) sep(isep)
            if (eof < 0) exit
          end do
          exit
        enddo
        close(172)
      else
        allocate (sep(0:0))
      end if

      db_mx%septic = imax

      return
      end subroutine sep_read
