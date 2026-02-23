       subroutine ch_read_parms_cal
      
       use calibration_data_module
       use input_file_module
       use utils, only: open_file

      implicit none

      character (len=80) :: titldum = "" !             |title of file
      character (len=80) :: header = ""  !             |header of file
      integer :: eof = 0               !             |end of file
      integer :: mchp = 0              !             |ending of loop
      integer :: i = 0                 !none         |counter

       eof = 0

      if (open_file(107, in_chg%ch_sed_parms_sft)) then
        do
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) mchp
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
          allocate (ch_prms(mchp))
          exit
        enddo
       
        do i = 1, mchp
          read (107,*,iostat=eof) ch_prms(i)%name, ch_prms(i)%chg_typ, ch_prms(i)%neg, ch_prms(i)%pos, ch_prms(i)%lo, ch_prms(i)%up
          if (eof < 0) exit 
        end do 
    
      else
           allocate (ch_prms(0:0))
      endif
      
      close(107)
      return
      end subroutine ch_read_parms_cal