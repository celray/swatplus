      subroutine gwflow_aqu_read

      use input_file_module
      use utils, only: open_file, to_lower

      implicit none

      character(len=512) :: line = ""
      character(len=25)  :: keyword = ""
      integer :: eof = 0
      integer :: ios = 0

      !! open gwflow.aqu (listed on aquifer line of file.cio)
      if (.not. open_file(107, in_aqu%gwflow)) return

      !! skip header line
      read (107, '(A)', iostat=eof) line
      if (eof < 0) then
        close (107)
        return
      end if

      !! read section lines: keyword  file1  file2  ...
      do
        read (107, '(A)', iostat=eof) line
        if (eof < 0) exit
        if (len_trim(line) == 0) cycle

        !! extract section keyword (first token)
        read (line, *, iostat=ios) keyword
        if (ios /= 0) cycle
        keyword = to_lower(trim(keyword))

        select case (trim(keyword))

        case ("basic")
          read (line, *, iostat=ios) keyword, in_gwf%gw_input

        case ("cells")
          read (line, *, iostat=ios) keyword, in_gwf%hrucell,       &
               in_gwf%lsucell, in_gwf%cellhru, in_gwf%huc12cell,   &
               in_gwf%con, in_gwf%chancells, in_gwf%rescells

        case ("exchange")
          read (line, *, iostat=ios) keyword, in_gwf%wetland,       &
               in_gwf%floodplain, in_gwf%canals

        case ("pumping")
          read (line, *, iostat=ios) keyword, in_gwf%pumpex,        &
               in_gwf%tiles

        case ("solutes")
          read (line, *, iostat=ios) keyword, in_gwf%solutes,       &
               in_gwf%solutes_minerals, in_gwf%streamobs

        case ("observations")
          read (line, *, iostat=ios) keyword, in_gwf%hru_pump_observe, &
               in_gwf%usgs_head

        end select
      end do

      close (107)
      return

      end subroutine gwflow_aqu_read
