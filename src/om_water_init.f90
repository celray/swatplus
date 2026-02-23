      subroutine om_water_init

      use basin_module
      use input_file_module
      use maximum_data_module
      use channel_data_module
      use hydrograph_module
      use sd_channel_module
      use constituent_mass_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""  !          |title of file
      character (len=80) :: header = "" !          |header of file
      integer :: eof = 0                !          |end of file
      integer :: imax = 0               !units     |description
      integer :: ichi = 0               !none      |counter

      eof = 0
      imax = 0

      if (open_file(105, in_init%om_water)) then
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

      db_mx%om_water_init = imax

      allocate (om_init_water(0:imax))
      allocate (om_init_name(0:imax))
      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit

       do ichi = 1, db_mx%om_water_init
         read (105,*,iostat=eof) titldum
         if (eof < 0) exit
         backspace (105)
         read (105,*,iostat=eof) om_init_name(ichi), om_init_water(ichi)
         if (eof < 0) exit
       end do
       close (105)
      exit
      enddo
      else
        allocate (om_init_water(0:0))
        allocate (om_init_name(0:0))
      endif

      return
      end subroutine om_water_init
