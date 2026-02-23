       subroutine aqu_read 
      
       use input_file_module
       use aquifer_module
       use basin_module !rtb gwflow
       use maximum_data_module
       use utils, only: open_file
       
       implicit none
      
       character (len=500) :: header = ""
       character (len=80) :: titldum = ""
       integer :: eof = 0         !                |end of file
       integer :: i = 0           !none            |counter
       integer :: imax = 0        !                |maximum count
       integer :: msh_aqp = 0     !none            |counter
       integer :: ish_aqp = 0     !none            |counter
       integer :: k = 0           !                |index
       
       msh_aqp = 0
       eof = 0
       imax = 0

       !! read shallow aquifer property data from aquifer.aqu
       if (open_file(107, in_aqu%aqu)) then
       do
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
            do while (eof == 0)
              read (107,*,iostat=eof) i
              if (eof < 0) exit
              imax = Max(imax,i)
              msh_aqp = msh_aqp + 1
            end do 
               
          db_mx%aqudb = msh_aqp
          allocate (aqudb(0:imax))
          rewind (107)
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
          
          do ish_aqp = 1, msh_aqp
            read (107,*,iostat=eof) i
            if (eof < 0) exit
            backspace (107)
            !! read from the aquifer database file named aquifer.aqu
            read (107,*,iostat=eof) k, aqudb(i)
            if (eof < 0) exit
          end do

          close (107)
          exit
          
          bsn_cc%gwflow = 0 ! rtb set gwflow module flag to 0
       enddo
       else
            allocate (aqudb(0:0))
       endif
          
       return
       end subroutine aqu_read         