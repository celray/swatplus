       subroutine ls_read_lsparms_cal

       use maximum_data_module
       use calibration_data_module
       use input_file_module
       use utils, only: open_file

       implicit none

       character (len=80) :: titldum = ""!           |title of file
       character (len=80) :: header = "" !           |header of file
       integer :: eof = 0              !           |end of file
       integer :: mlsp = 0             !none       |end of loop
       integer :: i = 0                !none       |counter

       eof = 0

      if (open_file(107, in_chg%wb_parms_sft)) then
       do
         read (107,*,iostat=eof) titldum
         if (eof < 0) exit
         read (107,*,iostat=eof) mlsp
         if (eof < 0) exit
         read (107,*,iostat=eof) header
         allocate (ls_prms(mlsp))
         if (eof < 0) exit
         exit
       enddo

       db_mx%lscal_prms = mlsp

       do i = 1, mlsp
         read (107,*,iostat=eof) ls_prms(i)%name, ls_prms(i)%chg_typ, ls_prms(i)%neg, ls_prms(i)%pos, ls_prms(i)%lo, ls_prms(i)%up
         if (eof < 0) exit
       end do
       close(107)
      else
        allocate (ls_prms(0:0))
      end if

      return
      end subroutine ls_read_lsparms_cal
