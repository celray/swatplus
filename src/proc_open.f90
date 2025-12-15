      subroutine proc_open
      
      use basin_module, only : pco
#ifdef SQLITE_ENABLED
      use sqlite_output_module, only : sqlite_init
#endif

      implicit none
      
      external :: header_aquifer, header_channel, header_const, header_hyd, header_lu_change, header_mgt, &
                  header_path, header_pest, header_reservoir, header_salt, header_sd_channel, header_snutc, &
                  header_water_allocation, header_wetland, header_write, header_yield, &
                  output_landscape_init, search

#ifdef SQLITE_ENABLED
      !! Initialize SQLite output if enabled
      if (pco%sqliteout == "y") then
        call sqlite_init("output.sqlite")
        call sqlite_create_output_tables()
      end if
#endif

      !! write headers in output files
      call output_landscape_init
      call header_channel
      call header_aquifer
      call header_sd_channel
      call header_mgt
      call header_lu_change
      call header_yield
      call header_hyd
      call header_reservoir
      call header_wetland
      call header_snutc
      call header_water_allocation
      
      call header_pest
      call header_path
      call header_salt !rtb salt
      call header_const !rtb cs

      call header_write
           
      return
      
      end subroutine proc_open
      
#ifdef SQLITE_ENABLED
      !> Create all output tables in SQLite database
      subroutine sqlite_create_output_tables()
        use sqlite_output_module
        
        implicit none
        
        character(len=16) :: col_base(7)
        character(len=16) :: typ_base(7)
        character(len=16) :: col_aqu(17)
        character(len=16) :: typ_aqu(17)
        
        !! Common columns for all output tables
        col_base = (/ "jday  ", "mon   ", "day   ", "yr    ", "gis_id", "unit  ", "name  " /)
        typ_base = (/ "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "INTEGER", "TEXT   " /)
        
        !! Aquifer balance columns
        col_aqu = (/ "flo     ", "dep_wt  ", "stor    ", "rchrg   ", "seep    ", "revap   ", &
                     "no3_rchg", "no3_loss", "no3_lat ", "no3_seep", "no3_st  ", "flo_cha ", &
                     "flo_res ", "flo_ls  ", "dep_wt2 ", "stor2   ", "no3_st2 " /)
        typ_aqu = (/ "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", &
                     "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", "REAL", "REAL" /)
        
        !! Create basin aquifer tables (daily, monthly, yearly, average annual)
        call sqlite_create_table("basin_aqu_day", [col_base, col_aqu], [typ_base, typ_aqu])
        call sqlite_create_table("basin_aqu_mon", [col_base, col_aqu], [typ_base, typ_aqu])
        call sqlite_create_table("basin_aqu_yr", [col_base, col_aqu], [typ_base, typ_aqu])
        call sqlite_create_table("basin_aqu_aa", [col_base, col_aqu], [typ_base, typ_aqu])
        
        !! Additional tables would be created here following the same pattern
        !! Each output type gets tables matching the text output file naming convention
        
      end subroutine sqlite_create_output_tables
#endif