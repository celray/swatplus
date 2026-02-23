      subroutine cs_plant_read !rtb cs
    
      use constituent_mass_module
      use input_file_module
      use maximum_data_module
      use cs_data_module
      use utils, only: open_file

      implicit none
 
      character (len=80) :: titldum = ""
      character (len=80) :: header = ""
      character (len=12) :: plant_name = ""
      integer :: iplant = 0


      !read plant boron tolerance data
      if (open_file(107, in_constit%plants_boron_cs)) then
        read(107,*) titldum

        !read plant boron tolerance parameters (a,b parameters from relative yield equations)
        read(107,*)
        read(107,*) bor_tol_sim !flag to simulate boron effect on plant growth 
        read(107,*) header
        read(107,*) header
        read(107,*) header
        read(107,*) header
        allocate (bor_stress_a(db_mx%plantparm), source = 0.)
        allocate (bor_stress_b(db_mx%plantparm), source = 0.)
        do iplant=1,db_mx%plantparm
          read(107,*) plant_name,bor_stress_a(iplant),bor_stress_b(iplant)
        enddo

        !close the file
        close(107)

      endif


      return
      end subroutine cs_plant_read