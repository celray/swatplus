      subroutine mgt_read_grazeops

      use input_file_module
      use maximum_data_module
      use mgt_operations_module
      use fertilizer_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: igrazop = 0          !none       |counter
      integer :: mgrazops = 0         !           |
      integer :: ifert = 0            !none       |counter

      mgrazops = 0
      eof = 0
      imax = 0

      !! read grazing operations
      if (open_file(107, in_ops%graze_ops)) then
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

        allocate (grazeop_db(0:imax))

        rewind (107)
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) header
        if (eof < 0) exit

        do igrazop = 1, imax
          read (107,*,iostat=eof) grazeop_db(igrazop)%name, grazeop_db(igrazop)%fertnm, grazeop_db(igrazop)%eat,    &
            grazeop_db(igrazop)%tramp, grazeop_db(igrazop)%manure, grazeop_db(igrazop)%biomin
          if (eof < 0) exit

          !xwalk fert name with fertilizer data base (fertilizer.frt)
          do ifert = 1, db_mx%fertparm
            if (grazeop_db(igrazop)%fertnm == fertdb(ifert)%fertnm) then
              grazeop_db(igrazop)%manure_id = ifert
              exit
            end if
          end do
        end do
      end do
      close(107)
      else
         allocate (grazeop_db(0:0))
      end if

      db_mx%grazeop_db = imax

      return
      end subroutine mgt_read_grazeops
