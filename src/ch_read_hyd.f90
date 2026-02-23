      subroutine ch_read_hyd
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use channel_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = "" !          |title of file
      character (len=80) :: header = ""  !          |header of file
      integer :: eof = 0               !          |end of file
      integer :: imax = 0              !units     |description
      integer :: ich = 0               !none      |counter

      eof = 0
      imax = 0

      if (open_file(105, in_cha%hyd)) then
      do
       read (105,*,iostat=eof) titldum
       if (eof < 0) exit
       read (105,*,iostat=eof) header
       if (eof < 0) exit
        do while (eof == 0)
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          imax = imax + 1
        end do
        
      db_mx%ch_hyd = imax
      
      allocate (ch_hyd(0:imax))
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
    
       do ich = 1, db_mx%ch_hyd
         read (105,*,iostat=eof) titldum
         if (eof < 0) exit
         backspace (105)
         read (105,*,iostat=eof) ch_hyd(ich)
         if (eof < 0) exit
         
        ch_hyd(ich)%alpha_bnk = Exp(-ch_hyd(ich)%alpha_bnk)
        if (ch_hyd(ich)%s <= 0.) ch_hyd(ich)%s = .0001
        if (ch_hyd(ich)%n <= 0.01) ch_hyd(ich)%n = .01
        if (ch_hyd(ich)%n >= 0.70) ch_hyd(ich)%n = 0.70
        if (ch_hyd(ich)%l <= 0.) ch_hyd(ich)%l = .0010
        if (ch_hyd(ich)%wdr <= 0.) ch_hyd(ich)%wdr = 3.5
        if (ch_hyd(ich)%side <= 1.e-6) ch_hyd(ich)%side = 2.0
        
       end do
       close (105)
      exit
      enddo
      else
        allocate (ch_hyd(0:0))
      endif

      return    
      end subroutine ch_read_hyd