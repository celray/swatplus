      subroutine cntbl_read

      use input_file_module
      use maximum_data_module
      use landuse_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !           |
      integer :: icno = 0             !none       |counter

      eof = 0
      imax = 0

    !! read all curve number data from cn.tbl
      if (open_file(107, in_lum%cntable_lum)) then
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

        allocate (cn(0:imax))

        rewind (107)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit

        do icno = 1, imax
          read (107,*,iostat=eof) cn(icno)
          if (eof < 0) exit
        end do
        exit
      enddo
      close(107)
      else
        allocate (cn(0:0))
      endif

      db_mx%cn_lu = imax

      return
      end subroutine cntbl_read
