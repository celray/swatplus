      subroutine ch_read_init
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use channel_data_module
      use sd_channel_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = "" !          |title of file
      character (len=80) :: header = ""  !          |header of file
      integer :: eof = 0               !          |end of file
      integer :: imax = 0              !units     |description
      integer :: ich = 0               !none      |counter

      eof = 0
      imax = 0

      if (open_file(105, in_cha%init)) then
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
        
      db_mx%ch_init = imax
      
      allocate (ch_init(0:imax))
      allocate (sd_init(0:imax))
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
      
       do ich = 1, db_mx%ch_init
         read (105,*,iostat=eof) ch_init(ich)
         if (eof < 0) exit
       end do
       close (105)
      exit
      enddo
      else
        allocate (ch_init(0:0))
        allocate (sd_init(0:0))
      endif

      return 
      end subroutine ch_read_init