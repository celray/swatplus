      subroutine cs_fert_read !rtb cs
       
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads constituent fertilizer loading (kg/ha) for various fertilizer types
    
      use constituent_mass_module
      use input_file_module
      use maximum_data_module
      use cs_module
      use utils, only: open_file

      implicit none
 
      character (len=80) :: titldum = ""
      character (len=80) :: header = ""
      !character (len=30) :: fert_name = ""
      integer :: icsi = 0
      integer :: eof = 0

      eof = 0
      
      !read constituent fertilizer loading (kg/ha)
      if (open_file(107, in_constit%fert_cs)) then
        read (107,*,iostat=eof) titldum
        read (107,*,iostat=eof) header
        
        !allocate fertilizer array
        allocate (fert_cs(db_mx%fertparm))
        
        !set flag
        fert_cs_flag = 1
        
        !read in the constituent fertilizer rates (kg/ha) for each fertilizer type
        do icsi=1,db_mx%fertparm
          read (107,*) fert_cs(icsi)  
        enddo
        close (107)
      end if
      
      return
      end subroutine cs_fert_read