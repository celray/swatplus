      subroutine res_read_csdb !rtb cs
      
      use input_file_module
      use maximum_data_module
      use reservoir_data_module
      use res_cs_module
      use utils, only: open_file

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads reservoir water quality parameters for constituents

      implicit none

      character (len=80) :: titldum = ""  !             |title of file
      character (len=80) :: header = "" !             |header of file
      integer :: i = 0                  !             |counter
      integer :: eof = 0                !             |end of file
      integer :: imax = 0               !             |determine max number for array (imax) and total number in file
      integer :: ires = 0               !none         |counter

      eof = 0
      imax = 0

      if (open_file(105, in_constit%res_cs)) then
      do
        read (105,*,iostat=eof) titldum
        read (105,*,iostat=eof) titldum
        do i=1,12
          read(105,*)
        enddo
        
        if (eof < 0) exit
        read (105,*,iostat=eof) header
        if (eof < 0) exit
          do while (eof == 0)
            read (105,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do   
          
        db_mx%res_cs = imax
        
        allocate (res_cs_data(0:imax))
        rewind (105)
        read (105,*,iostat=eof) titldum
        if (eof < 0) exit
        read (105,*,iostat=eof) titldum
        do i=1,12
          read(105,*)
        enddo
        read (105,*,iostat=eof) header
        if (eof < 0) exit
          
        do ires = 1, imax
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          backspace (105)
          read (105,*,iostat=eof) res_cs_data(ires)
          if (eof < 0) exit
        end do
        exit
      enddo
      else
        allocate (res_cs_data(0:0))
      endif
      close(105)
         
      return
      end subroutine res_read_csdb