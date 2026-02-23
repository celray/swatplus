      subroutine cal_parm_read

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this function computes new parameter value based on
!!    user defined change

      use input_file_module
      use maximum_data_module
      use calibration_data_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""                   !         |title of file
      character (len=80) :: header = ""                    !         |header of file
      integer :: eof = 0                                   !         |end of file
      integer :: imax = 0                                  !         |determine max number for array (imax) and total number in file
      integer :: mchg_par = 0                              !         |
      integer :: i = 0                                     !none     |counter

      imax = 0
      mchg_par = 0

      !!read parameter change values for calibration
      if (open_file(107, in_chg%cal_parms)) then
        do
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) mchg_par
          if (eof < 0) exit

          allocate (cal_parms(0:mchg_par))

          read (107,*,iostat=eof) header
          if (eof < 0) exit

          do i = 1, mchg_par
            read (107,*,iostat=eof) cal_parms(i)
            if (eof < 0) exit
          end do
          exit
        end do
        close (107)
      else
        allocate (cal_parms(0:0))
      end if

      db_mx%cal_parms = mchg_par

      return
      end subroutine cal_parm_read
