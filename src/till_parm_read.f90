      subroutine till_parm_read
      
      use input_file_module
      use maximum_data_module
      use tillage_data_module
      use utils, only: open_file
      
      implicit none      

      character (len=80) :: titldum = ""!           |title of file
      character (len=80) :: header = "" !           |header of file
      integer :: eof = 0              !           |end of file
      integer :: imax = 0             !none       |determine max number for array (imax) and total number in file
      integer :: itl = 0              !none       |counter
      integer :: mtl = 0              !           |
      
      eof = 0
      imax = 0
      mtl = 0
      bmix_idtill = 0
      
      if (open_file(105, in_parmdb%till_til)) then
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
            
          allocate (tilldb(0:imax))
          
          rewind (105)
          read (105,*,iostat=eof) titldum
          if (eof < 0) exit
          read (105,*,iostat=eof) header  
          if (eof < 0) exit
          
            do itl = 1, imax
              read (105,*,iostat=eof) tilldb(itl)
              if (tilldb(itl)%tillnm == "biomix") then
                bmix_idtill = itl
                bmix_eff = tilldb(itl)%effmix
                bmix_depth = tilldb(itl)%deptil
              endif
              if (eof < 0) exit
            end do    
          exit
        enddo
        if (bmix_idtill == 0) then
          bmix_eff = 0.2
          bmix_depth = 50.
        endif
      else
          allocate (tilldb(0:0))
      endif
      
      db_mx%tillparm = imax

      close (105)
      return
      end subroutine till_parm_read