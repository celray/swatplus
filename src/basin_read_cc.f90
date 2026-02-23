       subroutine basin_read_cc
      
       use input_file_module
       use basin_module
       use utils, only: open_file

       implicit none

       character (len=80) :: titldum = "" !             |title of file
       character (len=80) :: header = ""  !             |header
       integer :: eof = 0               !             |end of file

       eof = 0

       !! read basin
       if (open_file(107, in_basin%codes_bas)) then
       do
         read (107,*,iostat=eof) titldum
         if (eof < 0) exit
         read (107,*,iostat=eof) header
         if (eof < 0) exit
         read (107,*,iostat=eof) bsn_cc
         if (eof < 0) exit
         exit
       enddo
       close(107)
       endif

       if (bsn_cc%pet == 3) then
       if (open_file(140, in_cli%pet_cli)) then
       do
        read (140,*,iostat=eof) titldum
        if (eof < 0) exit
        read (140,*,iostat=eof) header
        if (eof < 0) exit
        read (140,*,iostat = eof) titldum
        exit
       end do
       end if
       end if
       
       return
      end subroutine basin_read_cc