       subroutine scen_read_filtstrip

       use input_file_module
       use maximum_data_module
       use mgt_operations_module
       use utils, only: open_file

       implicit none

       character (len=80) :: titldum = ""!           |title of file
       character (len=80) :: header = "" !           |header of file
       integer :: eof = 0              !           |end of file
       integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
       integer :: ifiltop = 0          !none       |counter

       eof = 0
       imax = 0

       !! read filter strip operations
       if (open_file(107, in_str%fstrip_str)) then
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

         allocate (filtstrip_db(0:imax))

         rewind (107)
         read (107,*,iostat=eof) titldum
         if (eof < 0) exit
         read (107,*,iostat=eof) header
         if (eof < 0) exit

         do ifiltop = 1, imax
           read (107,*,iostat=eof) filtstrip_db(ifiltop)
           if (eof < 0) exit
         end do
         exit
       enddo
       close(107)
       else
         allocate (filtstrip_db(0:0))
       endif

       db_mx%filtop_db = imax
       return
      end subroutine scen_read_filtstrip
