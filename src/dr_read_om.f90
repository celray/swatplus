      subroutine dr_read_om

      use dr_module
      use constituent_mass_module
      use hydrograph_module
      use input_file_module
      use organic_mineral_mass_module
      use maximum_data_module
      use utils, only: open_file

      implicit none
      
      character (len=80) :: titldum = ""
      character (len=80) :: header = ""
      integer :: eof = 0
      integer :: imax = 0
      integer :: ob1 = 0
      integer :: ob2 = 0
      integer :: ii = 0
      integer :: idr = 0
      integer :: idr_om = 0
      integer :: iob = 0

      eof = 0
      imax = 0
      
      !read delivery ratio organic-mineral data
      if (open_file(107, in_delr%om)) then
        do
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
          imax = 0
          do while (eof == 0)
            read (107,*,iostat=eof) titldum
            if (eof < 0) exit
            imax = imax + 1
          end do
          
          db_mx%dr_om = imax
          
          allocate (dr(0:imax))       !
          allocate (dr_om_num(0:imax), source = 0)
          allocate (dr_om_name(0:imax))
          rewind (107)
          read (107,*,iostat=eof) titldum
          if (eof < 0) exit
          read (107,*,iostat=eof) header
          if (eof < 0) exit
      
          !read all delivery ratio data
          do ii = 1, db_mx%dr_om
            read (107,*,iostat=eof) titldum
            if (eof < 0) exit
            backspace (107)
            read (107,*,iostat=eof) dr_om_name(ii), dr(ii)   
            if (eof < 0) exit
          end do
          close (107)
          exit
        end do
      end if
      
      ! xwalk with dr file to get sequential number
      do idr = 1, db_mx%dr
        do idr_om = 1, db_mx%dr_om
          if (dr_db(idr)%om_file == dr_om_name(idr_om)) then
            dr_om_num(idr) = idr_om
            exit
          end if
        end do
      end do
      
      !set dr_om object hydrograph
      ob1 = sp_ob1%dr
      ob2 = sp_ob1%dr + sp_ob%dr - 1
      do iob = ob1, ob2
        idr = ob(iob)%props
        idr_om = dr_om_num(idr)
        ob(iob)%hd(1) = dr(idr_om)
      end do
      
      return

      end subroutine dr_read_om