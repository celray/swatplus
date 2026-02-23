      subroutine mgt_read_mgtops

      use input_file_module
      use maximum_data_module
      use mgt_operations_module
      use utils, only: open_file

      implicit none

      external :: read_mgtops, mgt_operatn

      integer :: nops = 0             !           |end of loop
      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: iops = 0             !none       |counter
      integer :: nauto = 0            !           |end of loop
      integer :: iauto = 0            !none       |counter
      integer :: isched = 0           !none       |counter
      integer :: m_autos = 0          !           |end of loop

      eof = 0
      imax = 0

      !!   read mgtops.dat file
      !! calculate number of records in management
      if (open_file(107, in_lum%management_sch)) then
      do
       read (107,*,iostat=eof) titldum
       if (eof < 0) exit
       read (107,*,iostat=eof) header
       if (eof < 0) exit
       do while (eof == 0)
         read (107,*,iostat=eof) titldum, nops, nauto
         if (eof < 0) exit
         do iauto = 1, nauto
           read (107,*,iostat=eof) titldum
           if (eof < 0) exit
         end do
         do iops = 1, nops
           read (107,*,iostat=eof) titldum
           if (eof < 0) exit
         end do
         imax = imax + 1
       end do

       allocate (sched(0:imax))

       rewind (107)
       read (107,*,iostat=eof) titldum
       if (eof < 0) exit
       read (107,*,iostat=eof) header
       if (eof < 0) exit

       do isched = 1, imax
         read (107,*,iostat=eof)  sched(isched)%name, sched(isched)%num_ops, sched(isched)%num_autos
         if (eof < 0) exit
         !! allocate and read the auto operations
         m_autos = sched(isched)%num_autos
         if (m_autos > 0) then
           allocate (sched(isched)%auto_name(m_autos))
           allocate (sched(isched)%num_db(m_autos), source = 0)
           do iauto = 1, m_autos
             read (107,*,iostat=eof)  sched(isched)%auto_name(iauto)

             !! check to see if generic table - ie. plant-harv for single summer crop
             if (sched(isched)%auto_name(iauto) == "pl_hv_summer1" .or.      &
                 sched(isched)%auto_name(iauto) == "pl_hv_winter1") then
               allocate (sched(isched)%auto_crop(1))
               backspace (107)
               sched(isched)%auto_crop_num = 1
               read (107,*,iostat=eof)  sched(isched)%auto_name(iauto), sched(isched)%auto_crop
             end if
             if (sched(isched)%auto_name(iauto) == "pl_hv_summer2") then
               allocate (sched(isched)%auto_crop(2))
               backspace (107)
               sched(isched)%auto_crop_num = 1
               read (107,*,iostat=eof)  sched(isched)%auto_name(iauto), sched(isched)%auto_crop
             end if

             if (eof < 0) exit
           end do
         end if
         !! allocate and read the scheduled operations
         allocate (sched(isched)%mgt_ops(sched(isched)%num_ops))
         call read_mgtops(isched)
       end do
       exit
      enddo
      close(107)
      else
        allocate (sched(0:0))
      endif
      db_mx%mgt_ops = imax

      return
      end subroutine mgt_read_mgtops
