      !populate initial constituent data for aquifers
      subroutine aqu_read_init !rtb cs
      
      use basin_module
      use input_file_module
      use maximum_data_module
      use aquifer_module
      use aqu_pesticide_module
      use hydrograph_module
      use constituent_mass_module
      use utils, only: open_file
      
      implicit none      
      
      character (len=80) :: titldum = ""  !             |title of file
      character (len=80) :: header = "" !             |header of file
      integer :: eof = 0                !             |end of file
      integer :: imax = 0               !             |determine max number for array (imax) and total number in file
      integer :: iaqu = 0               !none         |counter
      integer :: iaq = 0
      integer :: ics = 0
      
      eof = 0
      imax = 0
      
      !read init
      if (open_file(105, in_aqu%init)) then
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

      allocate (aqu_init(0:imax))
      allocate (aqu_init_dat_c(0:imax))

      rewind (105)
      read (105,*,iostat=eof) titldum
      if (eof < 0) exit
      read (105,*,iostat=eof) header
      if (eof < 0) exit
           
       do iaqu = 1, imax
         read (105,*,iostat=eof) aqu_init_dat_c(iaqu)
         if (eof < 0) exit
       end do
       
       end do
       close (105)

      else
        allocate (aqu_init(0:0))
      end if

      !! initialize organics and constituents for each aquifer object
      do iaq = 1, sp_ob%aqu

        !! initial organic mineral
        do ics = 1, db_mx%om_water_init
          !! initializing organics in aqu_initial - do it here later
          !if (aqu_init(ini)%org_min == 0) write (9001,*) om_init_name(ics), " not found"
        end do
            
      end do

      return
      end subroutine aqu_read_init