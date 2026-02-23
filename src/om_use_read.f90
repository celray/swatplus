      subroutine om_use_read
      
      use input_file_module
      use water_allocation_module
      use mgt_operations_module
      use maximum_data_module
      use hydrograph_module
      use constituent_mass_module
      use utils, only: open_file

      implicit none

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      logical :: i_exist              !none       |check to determine if file exists
      integer :: iom_use = 0          !none       |

      eof = 0
      imax = 0

      !! read water allocation inputs

      if (open_file(107, in_wal%om_use)) then
      do
        read (107,*,iostat=eof) titldum
        if (eof < 0) exit
        read (107,*,iostat=eof) imax
        read (107,*,iostat=eof) header
        db_mx%om_use = imax
        if (eof < 0) exit
        
        allocate (wuse_om_efflu(imax))
        allocate (om_use_name(imax))

        do iom_use = 1, imax
          if (eof < 0) exit 
          read (107,*,iostat=eof) om_use_name(iom_use), wuse_om_efflu(iom_use)
        end do
      end do
      else
        allocate (wuse_om_efflu(0:0))
        allocate (om_use_name(0:0))
      end if

      close(107)

      return
      end subroutine om_use_read