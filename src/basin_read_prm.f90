      subroutine basin_read_prm
      
      use input_file_module
      use basin_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header
      integer :: eof = 0              !           |end of file

      eof = 0

      if (open_file(107, in_basin%parms_bas)) then
        !! read basin parameters
      do
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit
        read (107,*,iostat=eof) bsn_prm
        if (eof < 0) exit
        exit
      enddo
        close(107)
      end if
      
      return 
      end subroutine basin_read_prm