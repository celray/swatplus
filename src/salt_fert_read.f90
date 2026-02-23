      subroutine salt_fert_read !rtb salt
       
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine reads salt ion fertilizer loading (kg/ha) for various fertilizer types
  
      use constituent_mass_module
      use input_file_module
      use maximum_data_module
      use salt_module
      use utils, only: open_file

      implicit none
 
      character (len=80) :: titldum = ""
      character (len=80) :: header = ""
      integer :: isalti = 0
      integer :: eof = 0

      eof = 0
      
      !read salt fertilizer loading (kg/ha)
      if (open_file(107, in_salt%fert_slt)) then
        read (107,*,iostat=eof) titldum
        read (107,*,iostat=eof) header
        
        !allocate fertilizer array
        allocate (fert_salt(db_mx%fertparm))
        
        !set flag
        fert_salt_flag = 1
        
        !read in the salt ion fertilizer rates (kg/ha) for each fertilizer type
        do isalti=1,db_mx%fertparm
          read (107,*) fert_salt(isalti)  
        enddo
        close (107)
      end if
      
      return
      end subroutine salt_fert_read