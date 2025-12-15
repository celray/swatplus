#ifdef SQLITE_ENABLED
      module sqlite_output_module
      !! High-level SQLite output helper module for SWAT+
      !! Provides easy table creation and row insertion matching text output patterns
      
      use sqlite3_module
      use output_path_module, only: out_path
      
      implicit none
      
      private
      
      ! Public procedures - core functions
      public :: sqlite_init, sqlite_finalize
      public :: sqlite_create_table, sqlite_begin_transaction, sqlite_end_transaction
      public :: sqlite_insert_row_start, sqlite_insert_row_start_named
      public :: sqlite_insert_row_add_int, sqlite_insert_row_add_int8, sqlite_insert_row_add_real
      public :: sqlite_insert_row_add_text, sqlite_insert_row_execute
      public :: sqlite_db
      
      ! Public procedures - insert helpers for each output type
      public :: sqlite_insert_basin_aqu, sqlite_insert_aquifer
      public :: sqlite_insert_basin_wb, sqlite_insert_basin_nb, sqlite_insert_basin_ls, sqlite_insert_basin_pw
      public :: sqlite_insert_basin_ch, sqlite_insert_channel
      public :: sqlite_insert_hru_wb, sqlite_insert_hru_nb, sqlite_insert_hru_ls, sqlite_insert_hru_pw
      
      ! Global database handle
      type(sqlite3_db), save :: sqlite_db
      logical, save :: sqlite_initialized = .false.
      
      ! Current insert statement
      type(sqlite3_stmt_ptr), save :: current_stmt
      integer, save :: current_bind_idx = 0
      logical, save :: in_transaction = .false.
      
      ! Transaction batching
      integer, save :: rows_in_transaction = 0
      integer, parameter :: TRANSACTION_BATCH_SIZE = 10000
      
      ! Table creation tracking - stores names of tables already created
      integer, parameter :: MAX_TABLES = 200
      character(len=64), save :: created_tables(MAX_TABLES)
      integer, save :: num_created_tables = 0
      
      contains
      
      !> Check if a table has already been created
      logical function table_exists_in_registry(table_name)
        character(len=*), intent(in) :: table_name
        integer :: i
        
        table_exists_in_registry = .false.
        do i = 1, num_created_tables
          if (trim(created_tables(i)) == trim(table_name)) then
            table_exists_in_registry = .true.
            return
          end if
        end do
      end function table_exists_in_registry
      
      !> Add a table to the created tables registry
      subroutine register_table(table_name)
        character(len=*), intent(in) :: table_name
        
        if (num_created_tables < MAX_TABLES) then
          num_created_tables = num_created_tables + 1
          created_tables(num_created_tables) = trim(table_name)
        end if
      end subroutine register_table
      
      !> Initialize SQLite database
      subroutine sqlite_init(dbname)
        character(len=*), intent(in) :: dbname
        character(len=512) :: full_path
        integer :: rc
        
        if (sqlite_initialized) return
        
        ! Build full path with output directory
        if (len_trim(out_path) > 0) then
          full_path = trim(out_path) // trim(dbname)
        else
          full_path = trim(dbname)
        end if
        
        rc = sqlite3_open(trim(full_path), sqlite_db)
        if (rc /= SQLITE_OK) then
          write(*,*) "Error opening SQLite database: ", trim(full_path)
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
          return
        end if
        
        ! Enable WAL mode for better concurrent access
        rc = sqlite3_exec(sqlite_db, "PRAGMA journal_mode=WAL;")
        
        ! Optimize for bulk inserts
        rc = sqlite3_exec(sqlite_db, "PRAGMA synchronous=NORMAL;")
        rc = sqlite3_exec(sqlite_db, "PRAGMA cache_size=10000;")
        rc = sqlite3_exec(sqlite_db, "PRAGMA temp_store=MEMORY;")
        
        sqlite_initialized = .true.
        write(*,*) "SQLite output database initialized: ", trim(full_path)
        
      end subroutine sqlite_init
      
      !> Finalize SQLite database
      subroutine sqlite_finalize()
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        ! End any open transaction
        if (in_transaction) then
          rc = sqlite3_exec(sqlite_db, "COMMIT;")
          in_transaction = .false.
        end if
        
        rc = sqlite3_close(sqlite_db)
        sqlite_initialized = .false.
        
      end subroutine sqlite_finalize
      
      !> Create a table with columns from header string
      !! table_name: name of table (derived from output filename without extension)
      !! columns: array of column names
      !! types: array of column types ('INTEGER', 'REAL', 'TEXT')
      subroutine sqlite_create_table(table_name, columns, types)
        character(len=*), intent(in) :: table_name
        character(len=*), intent(in) :: columns(:)
        character(len=*), intent(in) :: types(:)
        character(len=4096) :: sql
        integer :: i, n, rc
        
        if (.not. sqlite_initialized) return
        
        n = size(columns)
        if (n /= size(types)) then
          write(*,*) "Error: columns and types arrays must have same size"
          return
        end if
        
        ! Build CREATE TABLE IF NOT EXISTS statement
        sql = "CREATE TABLE IF NOT EXISTS " // trim(table_name) // " ("
        
        do i = 1, n
          if (i > 1) sql = trim(sql) // ", "
          sql = trim(sql) // '"' // trim(adjustl(columns(i))) // '" ' // trim(types(i))
        end do
        
        sql = trim(sql) // ");"
        
        rc = sqlite3_exec(sqlite_db, trim(sql))
        if (rc /= SQLITE_OK) then
          write(*,*) "Error creating table: ", trim(table_name)
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
        end if
        
      end subroutine sqlite_create_table
      
      !> Begin a transaction (for bulk inserts)
      subroutine sqlite_begin_transaction()
        integer :: rc
        
        if (.not. sqlite_initialized) return
        if (in_transaction) return
        
        rc = sqlite3_exec(sqlite_db, "BEGIN TRANSACTION;")
        in_transaction = .true.
        rows_in_transaction = 0
        
      end subroutine sqlite_begin_transaction
      
      !> End a transaction
      subroutine sqlite_end_transaction()
        integer :: rc
        
        if (.not. sqlite_initialized) return
        if (.not. in_transaction) return
        
        rc = sqlite3_exec(sqlite_db, "COMMIT;")
        in_transaction = .false.
        rows_in_transaction = 0
        
      end subroutine sqlite_end_transaction
      
      !> Start building an INSERT statement with named columns
      !! Creates the table with proper column names if it doesn't exist
      subroutine sqlite_insert_row_start_named(table_name, col_names)
        character(len=*), intent(in) :: table_name
        character(len=*), intent(in) :: col_names(:)
        character(len=8192) :: sql
        integer :: i, n, rc
        
        if (.not. sqlite_initialized) return
        
        n = size(col_names)
        
        ! Auto-start transaction if not in one
        if (.not. in_transaction) then
          call sqlite_begin_transaction()
        end if
        
        ! Create table with proper column names if not already done
        if (.not. table_exists_in_registry(table_name)) then
          call sqlite_create_table_with_names(table_name, col_names)
          call register_table(table_name)
        end if
        
        ! Build parameterized INSERT statement
        sql = "INSERT INTO """ // trim(table_name) // """ VALUES (?"
        do i = 2, n
          sql = trim(sql) // ",?"
        end do
        sql = trim(sql) // ");"
        
        ! Prepare statement
        rc = sqlite3_prepare_v2(sqlite_db, trim(sql), current_stmt)
        if (rc /= SQLITE_OK) then
          write(*,*) "Error preparing insert for: ", trim(table_name)
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
        end if
        
        current_bind_idx = 0
        
      end subroutine sqlite_insert_row_start_named
      
      !> Create a table with proper column names
      subroutine sqlite_create_table_with_names(table_name, col_names)
        character(len=*), intent(in) :: table_name
        character(len=*), intent(in) :: col_names(:)
        character(len=8192) :: sql
        integer :: i, n, rc
        
        n = size(col_names)
        
        ! Build CREATE TABLE statement - use TEXT type for flexibility
        sql = "CREATE TABLE IF NOT EXISTS """ // trim(table_name) // """ ("""
        sql = trim(sql) // trim(adjustl(col_names(1))) // """ TEXT"
        do i = 2, n
          sql = trim(sql) // ", """ // trim(adjustl(col_names(i))) // """ TEXT"
        end do
        sql = trim(sql) // ");"
        
        rc = sqlite3_exec(sqlite_db, trim(sql))
        if (rc /= SQLITE_OK) then
          write(*,*) "Error creating table: ", trim(table_name)
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
        end if
        
      end subroutine sqlite_create_table_with_names
      
      !> Start building an INSERT statement (legacy - uses generic column names)
      !! Auto-creates the table if it doesn't exist
      subroutine sqlite_insert_row_start(table_name, num_cols)
        character(len=*), intent(in) :: table_name
        integer, intent(in) :: num_cols
        character(len=4096) :: sql
        integer :: i, rc
        
        if (.not. sqlite_initialized) return
        
        ! Auto-start transaction if not in one
        if (.not. in_transaction) then
          call sqlite_begin_transaction()
        end if
        
        ! Build parameterized INSERT statement
        sql = "INSERT INTO " // trim(table_name) // " VALUES (?"
        do i = 2, num_cols
          sql = trim(sql) // ",?"
        end do
        sql = trim(sql) // ");"
        
        ! Prepare statement
        rc = sqlite3_prepare_v2(sqlite_db, trim(sql), current_stmt)
        if (rc /= SQLITE_OK) then
          ! Table probably doesn't exist - create it with generic column names
          call sqlite_auto_create_table(table_name, num_cols)
          ! Retry prepare
          rc = sqlite3_prepare_v2(sqlite_db, trim(sql), current_stmt)
          if (rc /= SQLITE_OK) then
            write(*,*) "Error preparing insert for: ", trim(table_name)
            write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
          end if
        end if
        
        current_bind_idx = 0
        
      end subroutine sqlite_insert_row_start
      
      !> Auto-create a table with generic column names when first accessed
      subroutine sqlite_auto_create_table(table_name, num_cols)
        character(len=*), intent(in) :: table_name
        integer, intent(in) :: num_cols
        character(len=4096) :: sql
        character(len=8) :: col_name
        integer :: i, rc
        
        ! Build CREATE TABLE statement with generic column names
        ! Use BLOB type which accepts any SQLite type
        sql = "CREATE TABLE IF NOT EXISTS " // trim(table_name) // " (c1 BLOB"
        do i = 2, num_cols
          write(col_name, '("c",I0)') i
          sql = trim(sql) // ", " // trim(col_name) // " BLOB"
        end do
        sql = trim(sql) // ");"
        
        rc = sqlite3_exec(sqlite_db, trim(sql))
        if (rc /= SQLITE_OK) then
          write(*,*) "Error auto-creating table: ", trim(table_name)
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
        end if
        
      end subroutine sqlite_auto_create_table
      
      !> Add integer value to current insert
      subroutine sqlite_insert_row_add_int(val)
        integer, intent(in) :: val
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        current_bind_idx = current_bind_idx + 1
        rc = sqlite3_bind_int(current_stmt, current_bind_idx, val)
        
      end subroutine sqlite_insert_row_add_int
      
      !> Add 64-bit integer value to current insert
      subroutine sqlite_insert_row_add_int8(val)
        integer*8, intent(in) :: val
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        current_bind_idx = current_bind_idx + 1
        rc = sqlite3_bind_int64(current_stmt, current_bind_idx, val)
        
      end subroutine sqlite_insert_row_add_int8
      
      !> Add real value to current insert
      subroutine sqlite_insert_row_add_real(val)
        real, intent(in) :: val
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        current_bind_idx = current_bind_idx + 1
        rc = sqlite3_bind_double(current_stmt, current_bind_idx, real(val, 8))
        
      end subroutine sqlite_insert_row_add_real
      
      !> Add text value to current insert
      subroutine sqlite_insert_row_add_text(val)
        character(len=*), intent(in) :: val
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        current_bind_idx = current_bind_idx + 1
        rc = sqlite3_bind_text(current_stmt, current_bind_idx, trim(val))
        
      end subroutine sqlite_insert_row_add_text
      
      !> Execute the current insert statement
      subroutine sqlite_insert_row_execute()
        integer :: rc
        
        if (.not. sqlite_initialized) return
        
        rc = sqlite3_step(current_stmt)
        if (rc /= SQLITE_DONE .and. rc /= SQLITE_ROW) then
          write(*,*) "Error executing insert, rc=", rc
          write(*,*) "SQLite error: ", trim(sqlite3_errmsg(sqlite_db))
        end if
        
        rc = sqlite3_finalize(current_stmt)
        current_bind_idx = 0
        
        ! Auto-commit transaction after batch size
        rows_in_transaction = rows_in_transaction + 1
        if (rows_in_transaction >= TRANSACTION_BATCH_SIZE) then
          call sqlite_end_transaction()
          call sqlite_begin_transaction()
        end if
        
      end subroutine sqlite_insert_row_execute
      
      !> Insert basin aquifer data into SQLite
      !! Helper subroutine for basin_aquifer_output
      subroutine sqlite_insert_basin_aqu(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, aqu_data)
        use aquifer_module, only : aquifer_dynamic
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(aquifer_dynamic), intent(in) :: aqu_data
        character(len=32) :: table_name
        
        character(len=16) :: col_names(24)
        
        if (.not. sqlite_initialized) return
        
        ! Determine table name based on period (matches text output filename pattern)
        table_name = "basin_aqu_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "flo", "dep_wt", "stor", "rchrg", "seep", "revap", "no3_rchg", "no3_loss", &
                     "no3_lat", "no3_seep", "no3_st", "flo_cha", "flo_res", "flo_ls", "minp", "cbn", "orgn"]
        
        ! Insert row with all aquifer dynamic fields
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(aqu_data%flo)
        call sqlite_insert_row_add_real(aqu_data%dep_wt)
        call sqlite_insert_row_add_real(aqu_data%stor)
        call sqlite_insert_row_add_real(aqu_data%rchrg)
        call sqlite_insert_row_add_real(aqu_data%seep)
        call sqlite_insert_row_add_real(aqu_data%revap)
        call sqlite_insert_row_add_real(aqu_data%no3_rchg)
        call sqlite_insert_row_add_real(aqu_data%no3_loss)
        call sqlite_insert_row_add_real(aqu_data%no3_lat)
        call sqlite_insert_row_add_real(aqu_data%no3_seep)
        call sqlite_insert_row_add_real(aqu_data%no3_st)
        call sqlite_insert_row_add_real(aqu_data%flo_cha)
        call sqlite_insert_row_add_real(aqu_data%flo_res)
        call sqlite_insert_row_add_real(aqu_data%flo_ls)
        call sqlite_insert_row_add_real(aqu_data%minp)
        call sqlite_insert_row_add_real(aqu_data%cbn)
        call sqlite_insert_row_add_real(aqu_data%orgn)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_aqu
      
      !> Insert aquifer data into SQLite (for individual aquifer output)
      subroutine sqlite_insert_aquifer(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, aqu_data)
        use aquifer_module, only : aquifer_dynamic
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(aquifer_dynamic), intent(in) :: aqu_data
        character(len=32) :: table_name
        character(len=16) :: col_names(24)
        
        if (.not. sqlite_initialized) return
        
        table_name = "aquifer_" // trim(period)
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "flo", "dep_wt", "stor", "rchrg", "seep", "revap", "no3_rchg", "no3_loss", &
                     "no3_lat", "no3_seep", "no3_st", "flo_cha", "flo_res", "flo_ls", "minp", "cbn", "orgn"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(aqu_data%flo)
        call sqlite_insert_row_add_real(aqu_data%dep_wt)
        call sqlite_insert_row_add_real(aqu_data%stor)
        call sqlite_insert_row_add_real(aqu_data%rchrg)
        call sqlite_insert_row_add_real(aqu_data%seep)
        call sqlite_insert_row_add_real(aqu_data%revap)
        call sqlite_insert_row_add_real(aqu_data%no3_rchg)
        call sqlite_insert_row_add_real(aqu_data%no3_loss)
        call sqlite_insert_row_add_real(aqu_data%no3_lat)
        call sqlite_insert_row_add_real(aqu_data%no3_seep)
        call sqlite_insert_row_add_real(aqu_data%no3_st)
        call sqlite_insert_row_add_real(aqu_data%flo_cha)
        call sqlite_insert_row_add_real(aqu_data%flo_res)
        call sqlite_insert_row_add_real(aqu_data%flo_ls)
        call sqlite_insert_row_add_real(aqu_data%minp)
        call sqlite_insert_row_add_real(aqu_data%cbn)
        call sqlite_insert_row_add_real(aqu_data%orgn)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_aquifer
      
      !> Insert basin water balance data into SQLite
      subroutine sqlite_insert_basin_wb(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, wb_data)
        use output_landscape_module, only : output_waterbal
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(output_waterbal), intent(in) :: wb_data
        character(len=32) :: table_name
        character(len=16) :: col_names(53)
        
        if (.not. sqlite_initialized) return
        
        table_name = "basin_wb_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "precip", "snofall", "snomlt", "surq_gen", "latq", "wateryld", "perc", "et", &
                     "ecanopy", "eplant", "esoil", "surq_cont", "cn", "sw_init", "sw_final", "sw", &
                     "sw_300", "sno_init", "sno_final", "snopack", "pet", "qtile", "irr", "surq_runon", &
                     "latq_runon", "overbank", "surq_cha", "surq_res", "surq_ls", "latq_cha", "latq_res", &
                     "latq_ls", "gwsoil", "satex", "satex_chan", "delsw", "lagsurf", "laglatq", "lagsatex", &
                     "wet_evap", "wet_out", "wet_stor", "null1", "null2", "null3", "null4"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(wb_data%precip)
        call sqlite_insert_row_add_real(wb_data%snofall)
        call sqlite_insert_row_add_real(wb_data%snomlt)
        call sqlite_insert_row_add_real(wb_data%surq_gen)
        call sqlite_insert_row_add_real(wb_data%latq)
        call sqlite_insert_row_add_real(wb_data%wateryld)
        call sqlite_insert_row_add_real(wb_data%perc)
        call sqlite_insert_row_add_real(wb_data%et)
        call sqlite_insert_row_add_real(wb_data%ecanopy)
        call sqlite_insert_row_add_real(wb_data%eplant)
        call sqlite_insert_row_add_real(wb_data%esoil)
        call sqlite_insert_row_add_real(wb_data%surq_cont)
        call sqlite_insert_row_add_real(wb_data%cn)
        call sqlite_insert_row_add_real(wb_data%sw_init)
        call sqlite_insert_row_add_real(wb_data%sw_final)
        call sqlite_insert_row_add_real(wb_data%sw)
        call sqlite_insert_row_add_real(wb_data%sw_300)
        call sqlite_insert_row_add_real(wb_data%sno_init)
        call sqlite_insert_row_add_real(wb_data%sno_final)
        call sqlite_insert_row_add_real(wb_data%snopack)
        call sqlite_insert_row_add_real(wb_data%pet)
        call sqlite_insert_row_add_real(wb_data%qtile)
        call sqlite_insert_row_add_real(wb_data%irr)
        call sqlite_insert_row_add_real(wb_data%surq_runon)
        call sqlite_insert_row_add_real(wb_data%latq_runon)
        call sqlite_insert_row_add_real(wb_data%overbank)
        call sqlite_insert_row_add_real(wb_data%surq_cha)
        call sqlite_insert_row_add_real(wb_data%surq_res)
        call sqlite_insert_row_add_real(wb_data%surq_ls)
        call sqlite_insert_row_add_real(wb_data%latq_cha)
        call sqlite_insert_row_add_real(wb_data%latq_res)
        call sqlite_insert_row_add_real(wb_data%latq_ls)
        call sqlite_insert_row_add_real(wb_data%gwsoil)
        call sqlite_insert_row_add_real(wb_data%satex)
        call sqlite_insert_row_add_real(wb_data%satex_chan)
        call sqlite_insert_row_add_real(wb_data%delsw)
        call sqlite_insert_row_add_real(wb_data%lagsurf)
        call sqlite_insert_row_add_real(wb_data%laglatq)
        call sqlite_insert_row_add_real(wb_data%lagsatex)
        call sqlite_insert_row_add_real(wb_data%wet_evap)
        call sqlite_insert_row_add_real(wb_data%wet_out)
        call sqlite_insert_row_add_real(wb_data%wet_stor)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_wb
      
      !> Insert basin nutrient balance data into SQLite
      subroutine sqlite_insert_basin_nb(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, nb_data)
        use output_landscape_module, only : output_nutbal
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(output_nutbal), intent(in) :: nb_data
        character(len=32) :: table_name
        character(len=16) :: col_names(27)
        
        if (.not. sqlite_initialized) return
        
        table_name = "basin_nb_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "grazn", "grazp", "lab_min_p", "act_sta_p", "fertn", "fertp", "fixn", "denit", &
                     "act_nit_n", "act_sta_n", "org_lab_p", "rsd_nitorg_n", "rsd_laborg_p", "no3atmo", &
                     "nh4atmo", "nuptake", "puptake", "gwsoiln", "gwsoilp", "null1"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(nb_data%grazn)
        call sqlite_insert_row_add_real(nb_data%grazp)
        call sqlite_insert_row_add_real(nb_data%lab_min_p)
        call sqlite_insert_row_add_real(nb_data%act_sta_p)
        call sqlite_insert_row_add_real(nb_data%fertn)
        call sqlite_insert_row_add_real(nb_data%fertp)
        call sqlite_insert_row_add_real(nb_data%fixn)
        call sqlite_insert_row_add_real(nb_data%denit)
        call sqlite_insert_row_add_real(nb_data%act_nit_n)
        call sqlite_insert_row_add_real(nb_data%act_sta_n)
        call sqlite_insert_row_add_real(nb_data%org_lab_p)
        call sqlite_insert_row_add_real(nb_data%rsd_nitorg_n)
        call sqlite_insert_row_add_real(nb_data%rsd_laborg_p)
        call sqlite_insert_row_add_real(nb_data%no3atmo)
        call sqlite_insert_row_add_real(nb_data%nh4atmo)
        call sqlite_insert_row_add_real(nb_data%nuptake)
        call sqlite_insert_row_add_real(nb_data%puptake)
        call sqlite_insert_row_add_real(nb_data%gwsoiln)
        call sqlite_insert_row_add_real(nb_data%gwsoilp)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_nb
      
      !> Insert basin losses data into SQLite
      subroutine sqlite_insert_basin_ls(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, ls_data)
        use output_landscape_module, only : output_losses
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(output_losses), intent(in) :: ls_data
        character(len=32) :: table_name
        character(len=16) :: col_names(19)
        
        if (.not. sqlite_initialized) return
        
        table_name = "basin_ls_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "sedyld", "sedorgn", "sedorgp", "surqno3", "latno3", "surqsolp", "usle", &
                     "sedminp", "tileno3", "lchlabp", "tilelabp", "satexn"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(ls_data%sedyld)
        call sqlite_insert_row_add_real(ls_data%sedorgn)
        call sqlite_insert_row_add_real(ls_data%sedorgp)
        call sqlite_insert_row_add_real(ls_data%surqno3)
        call sqlite_insert_row_add_real(ls_data%latno3)
        call sqlite_insert_row_add_real(ls_data%surqsolp)
        call sqlite_insert_row_add_real(ls_data%usle)
        call sqlite_insert_row_add_real(ls_data%sedminp)
        call sqlite_insert_row_add_real(ls_data%tileno3)
        call sqlite_insert_row_add_real(ls_data%lchlabp)
        call sqlite_insert_row_add_real(ls_data%tilelabp)
        call sqlite_insert_row_add_real(ls_data%satexn)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_ls
      
      !> Insert basin plant weather data into SQLite
      subroutine sqlite_insert_basin_pw(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, pw_data)
        use output_landscape_module, only : output_plantweather
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(output_plantweather), intent(in) :: pw_data
        character(len=32) :: table_name
        character(len=16) :: col_names(32)
        
        if (.not. sqlite_initialized) return
        
        table_name = "basin_pw_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "lai", "bioms", "yield", "residue", "sol_tmp", "strsw", "strsa", "strstmp", &
                     "strsn", "strsp", "strss", "nplnt", "percn", "pplnt", "tmx", "tmn", "tmpav", &
                     "solrad", "wndspd", "rhum", "phubase0", "lai_max", "bm_max", "bm_grow", "c_gro"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(pw_data%lai)
        call sqlite_insert_row_add_real(pw_data%bioms)
        call sqlite_insert_row_add_real(pw_data%yield)
        call sqlite_insert_row_add_real(pw_data%residue)
        call sqlite_insert_row_add_real(pw_data%sol_tmp)
        call sqlite_insert_row_add_real(pw_data%strsw)
        call sqlite_insert_row_add_real(pw_data%strsa)
        call sqlite_insert_row_add_real(pw_data%strstmp)
        call sqlite_insert_row_add_real(pw_data%strsn)
        call sqlite_insert_row_add_real(pw_data%strsp)
        call sqlite_insert_row_add_real(pw_data%strss)
        call sqlite_insert_row_add_real(pw_data%nplnt)
        call sqlite_insert_row_add_real(pw_data%percn)
        call sqlite_insert_row_add_real(pw_data%pplnt)
        call sqlite_insert_row_add_real(pw_data%tmx)
        call sqlite_insert_row_add_real(pw_data%tmn)
        call sqlite_insert_row_add_real(pw_data%tmpav)
        call sqlite_insert_row_add_real(pw_data%solrad)
        call sqlite_insert_row_add_real(pw_data%wndspd)
        call sqlite_insert_row_add_real(pw_data%rhum)
        call sqlite_insert_row_add_real(pw_data%phubase0)
        call sqlite_insert_row_add_real(pw_data%lai_max)
        call sqlite_insert_row_add_real(pw_data%bm_max)
        call sqlite_insert_row_add_real(pw_data%bm_grow)
        call sqlite_insert_row_add_real(pw_data%c_gro)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_pw
      
      !> Insert basin channel data into SQLite
      subroutine sqlite_insert_basin_ch(period, jday, mon, day_mo, yrc, gis_id, unit_id, name, ch_data)
        use channel_module, only : ch_output
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(ch_output), intent(in) :: ch_data
        character(len=32) :: table_name
        character(len=16) :: col_names(61)
        
        if (.not. sqlite_initialized) return
        
        table_name = "basin_cha_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "gis_id", "unit", "name", &
                     "flo_in", "flo_out", "evap", "tloss", "sed_in", "sed_out", "sed_conc", &
                     "orgn_in", "orgn_out", "orgp_in", "orgp_out", "no3_in", "no3_out", "nh4_in", &
                     "nh4_out", "no2_in", "no2_out", "solp_in", "solp_out", "chla_in", "chla_out", &
                     "cbod_in", "cbod_out", "dis_in", "dis_out", "solpst_in", "solpst_out", &
                     "sorbpst_in", "sorbpst_out", "react", "volat", "setlpst", "resuspst", "difus", &
                     "reactb", "bury", "sedpest", "bacp", "baclp", "met1", "met2", "met3", &
                     "sand_in", "sand_out", "silt_in", "silt_out", "clay_in", "clay_out", &
                     "smag_in", "smag_out", "lag_in", "lag_out", "grvl_in", "grvl_out"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(ch_data%flo_in)
        call sqlite_insert_row_add_real(ch_data%flo_out)
        call sqlite_insert_row_add_real(ch_data%evap)
        call sqlite_insert_row_add_real(ch_data%tloss)
        call sqlite_insert_row_add_real(ch_data%sed_in)
        call sqlite_insert_row_add_real(ch_data%sed_out)
        call sqlite_insert_row_add_real(ch_data%sed_conc)
        call sqlite_insert_row_add_real(ch_data%orgn_in)
        call sqlite_insert_row_add_real(ch_data%orgn_out)
        call sqlite_insert_row_add_real(ch_data%orgp_in)
        call sqlite_insert_row_add_real(ch_data%orgp_out)
        call sqlite_insert_row_add_real(ch_data%no3_in)
        call sqlite_insert_row_add_real(ch_data%no3_out)
        call sqlite_insert_row_add_real(ch_data%nh4_in)
        call sqlite_insert_row_add_real(ch_data%nh4_out)
        call sqlite_insert_row_add_real(ch_data%no2_in)
        call sqlite_insert_row_add_real(ch_data%no2_out)
        call sqlite_insert_row_add_real(ch_data%solp_in)
        call sqlite_insert_row_add_real(ch_data%solp_out)
        call sqlite_insert_row_add_real(ch_data%chla_in)
        call sqlite_insert_row_add_real(ch_data%chla_out)
        call sqlite_insert_row_add_real(ch_data%cbod_in)
        call sqlite_insert_row_add_real(ch_data%cbod_out)
        call sqlite_insert_row_add_real(ch_data%dis_in)
        call sqlite_insert_row_add_real(ch_data%dis_out)
        call sqlite_insert_row_add_real(ch_data%solpst_in)
        call sqlite_insert_row_add_real(ch_data%solpst_out)
        call sqlite_insert_row_add_real(ch_data%sorbpst_in)
        call sqlite_insert_row_add_real(ch_data%sorbpst_out)
        call sqlite_insert_row_add_real(ch_data%react)
        call sqlite_insert_row_add_real(ch_data%volat)
        call sqlite_insert_row_add_real(ch_data%setlpst)
        call sqlite_insert_row_add_real(ch_data%resuspst)
        call sqlite_insert_row_add_real(ch_data%difus)
        call sqlite_insert_row_add_real(ch_data%reactb)
        call sqlite_insert_row_add_real(ch_data%bury)
        call sqlite_insert_row_add_real(ch_data%sedpest)
        call sqlite_insert_row_add_real(ch_data%bacp)
        call sqlite_insert_row_add_real(ch_data%baclp)
        call sqlite_insert_row_add_real(ch_data%met1)
        call sqlite_insert_row_add_real(ch_data%met2)
        call sqlite_insert_row_add_real(ch_data%met3)
        call sqlite_insert_row_add_real(ch_data%sand_in)
        call sqlite_insert_row_add_real(ch_data%sand_out)
        call sqlite_insert_row_add_real(ch_data%silt_in)
        call sqlite_insert_row_add_real(ch_data%silt_out)
        call sqlite_insert_row_add_real(ch_data%clay_in)
        call sqlite_insert_row_add_real(ch_data%clay_out)
        call sqlite_insert_row_add_real(ch_data%smag_in)
        call sqlite_insert_row_add_real(ch_data%smag_out)
        call sqlite_insert_row_add_real(ch_data%lag_in)
        call sqlite_insert_row_add_real(ch_data%lag_out)
        call sqlite_insert_row_add_real(ch_data%grvl_in)
        call sqlite_insert_row_add_real(ch_data%grvl_out)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_basin_ch
      
      !> Insert channel data into SQLite (individual channel)
      subroutine sqlite_insert_channel(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, ch_data)
        use channel_module, only : ch_output
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name
        type(ch_output), intent(in) :: ch_data
        character(len=32) :: table_name
        character(len=16) :: col_names(61)
        
        if (.not. sqlite_initialized) return
        
        table_name = "channel_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "flo_in", "flo_out", "evap", "tloss", "sed_in", "sed_out", "sed_conc", &
                     "orgn_in", "orgn_out", "orgp_in", "orgp_out", "no3_in", "no3_out", "nh4_in", &
                     "nh4_out", "no2_in", "no2_out", "solp_in", "solp_out", "chla_in", "chla_out", &
                     "cbod_in", "cbod_out", "dis_in", "dis_out", "solpst_in", "solpst_out", &
                     "sorbpst_in", "sorbpst_out", "react", "volat", "setlpst", "resuspst", "difus", &
                     "reactb", "bury", "sedpest", "bacp", "baclp", "met1", "met2", "met3", &
                     "sand_in", "sand_out", "silt_in", "silt_out", "clay_in", "clay_out", &
                     "smag_in", "smag_out", "lag_in", "lag_out", "grvl_in", "grvl_out"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(ch_data%flo_in)
        call sqlite_insert_row_add_real(ch_data%flo_out)
        call sqlite_insert_row_add_real(ch_data%evap)
        call sqlite_insert_row_add_real(ch_data%tloss)
        call sqlite_insert_row_add_real(ch_data%sed_in)
        call sqlite_insert_row_add_real(ch_data%sed_out)
        call sqlite_insert_row_add_real(ch_data%sed_conc)
        call sqlite_insert_row_add_real(ch_data%orgn_in)
        call sqlite_insert_row_add_real(ch_data%orgn_out)
        call sqlite_insert_row_add_real(ch_data%orgp_in)
        call sqlite_insert_row_add_real(ch_data%orgp_out)
        call sqlite_insert_row_add_real(ch_data%no3_in)
        call sqlite_insert_row_add_real(ch_data%no3_out)
        call sqlite_insert_row_add_real(ch_data%nh4_in)
        call sqlite_insert_row_add_real(ch_data%nh4_out)
        call sqlite_insert_row_add_real(ch_data%no2_in)
        call sqlite_insert_row_add_real(ch_data%no2_out)
        call sqlite_insert_row_add_real(ch_data%solp_in)
        call sqlite_insert_row_add_real(ch_data%solp_out)
        call sqlite_insert_row_add_real(ch_data%chla_in)
        call sqlite_insert_row_add_real(ch_data%chla_out)
        call sqlite_insert_row_add_real(ch_data%cbod_in)
        call sqlite_insert_row_add_real(ch_data%cbod_out)
        call sqlite_insert_row_add_real(ch_data%dis_in)
        call sqlite_insert_row_add_real(ch_data%dis_out)
        call sqlite_insert_row_add_real(ch_data%solpst_in)
        call sqlite_insert_row_add_real(ch_data%solpst_out)
        call sqlite_insert_row_add_real(ch_data%sorbpst_in)
        call sqlite_insert_row_add_real(ch_data%sorbpst_out)
        call sqlite_insert_row_add_real(ch_data%react)
        call sqlite_insert_row_add_real(ch_data%volat)
        call sqlite_insert_row_add_real(ch_data%setlpst)
        call sqlite_insert_row_add_real(ch_data%resuspst)
        call sqlite_insert_row_add_real(ch_data%difus)
        call sqlite_insert_row_add_real(ch_data%reactb)
        call sqlite_insert_row_add_real(ch_data%bury)
        call sqlite_insert_row_add_real(ch_data%sedpest)
        call sqlite_insert_row_add_real(ch_data%bacp)
        call sqlite_insert_row_add_real(ch_data%baclp)
        call sqlite_insert_row_add_real(ch_data%met1)
        call sqlite_insert_row_add_real(ch_data%met2)
        call sqlite_insert_row_add_real(ch_data%met3)
        call sqlite_insert_row_add_real(ch_data%sand_in)
        call sqlite_insert_row_add_real(ch_data%sand_out)
        call sqlite_insert_row_add_real(ch_data%silt_in)
        call sqlite_insert_row_add_real(ch_data%silt_out)
        call sqlite_insert_row_add_real(ch_data%clay_in)
        call sqlite_insert_row_add_real(ch_data%clay_out)
        call sqlite_insert_row_add_real(ch_data%smag_in)
        call sqlite_insert_row_add_real(ch_data%smag_out)
        call sqlite_insert_row_add_real(ch_data%lag_in)
        call sqlite_insert_row_add_real(ch_data%lag_out)
        call sqlite_insert_row_add_real(ch_data%grvl_in)
        call sqlite_insert_row_add_real(ch_data%grvl_out)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_channel
      
      !> Insert HRU water balance data into SQLite
      subroutine sqlite_insert_hru_wb(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, wb_data, plant_cov, mgt_ops)
        use output_landscape_module, only : output_waterbal
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name, plant_cov, mgt_ops
        type(output_waterbal), intent(in) :: wb_data
        character(len=32) :: table_name
        character(len=16) :: col_names(55)
        
        if (.not. sqlite_initialized) return
        
        table_name = "hru_wb_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "precip", "snofall", "snomlt", "surq_gen", "latq", "wateryld", "perc", "et", &
                     "ecanopy", "eplant", "esoil", "surq_cont", "cn", "sw_init", "sw_final", "sw", &
                     "sw_300", "sno_init", "sno_final", "snopack", "pet", "qtile", "irr", "surq_runon", &
                     "latq_runon", "overbank", "surq_cha", "surq_res", "surq_ls", "latq_cha", "latq_res", &
                     "latq_ls", "gwsoil", "satex", "satex_chan", "delsw", "lagsurf", "laglatq", "lagsatex", &
                     "wet_evap", "wet_out", "wet_stor", "plant_cov", "mgt_ops", "null1", "null2", "null3", "null4"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(wb_data%precip)
        call sqlite_insert_row_add_real(wb_data%snofall)
        call sqlite_insert_row_add_real(wb_data%snomlt)
        call sqlite_insert_row_add_real(wb_data%surq_gen)
        call sqlite_insert_row_add_real(wb_data%latq)
        call sqlite_insert_row_add_real(wb_data%wateryld)
        call sqlite_insert_row_add_real(wb_data%perc)
        call sqlite_insert_row_add_real(wb_data%et)
        call sqlite_insert_row_add_real(wb_data%ecanopy)
        call sqlite_insert_row_add_real(wb_data%eplant)
        call sqlite_insert_row_add_real(wb_data%esoil)
        call sqlite_insert_row_add_real(wb_data%surq_cont)
        call sqlite_insert_row_add_real(wb_data%cn)
        call sqlite_insert_row_add_real(wb_data%sw_init)
        call sqlite_insert_row_add_real(wb_data%sw_final)
        call sqlite_insert_row_add_real(wb_data%sw)
        call sqlite_insert_row_add_real(wb_data%sw_300)
        call sqlite_insert_row_add_real(wb_data%sno_init)
        call sqlite_insert_row_add_real(wb_data%sno_final)
        call sqlite_insert_row_add_real(wb_data%snopack)
        call sqlite_insert_row_add_real(wb_data%pet)
        call sqlite_insert_row_add_real(wb_data%qtile)
        call sqlite_insert_row_add_real(wb_data%irr)
        call sqlite_insert_row_add_real(wb_data%surq_runon)
        call sqlite_insert_row_add_real(wb_data%latq_runon)
        call sqlite_insert_row_add_real(wb_data%overbank)
        call sqlite_insert_row_add_real(wb_data%surq_cha)
        call sqlite_insert_row_add_real(wb_data%surq_res)
        call sqlite_insert_row_add_real(wb_data%surq_ls)
        call sqlite_insert_row_add_real(wb_data%latq_cha)
        call sqlite_insert_row_add_real(wb_data%latq_res)
        call sqlite_insert_row_add_real(wb_data%latq_ls)
        call sqlite_insert_row_add_real(wb_data%gwsoil)
        call sqlite_insert_row_add_real(wb_data%satex)
        call sqlite_insert_row_add_real(wb_data%satex_chan)
        call sqlite_insert_row_add_real(wb_data%delsw)
        call sqlite_insert_row_add_real(wb_data%lagsurf)
        call sqlite_insert_row_add_real(wb_data%laglatq)
        call sqlite_insert_row_add_real(wb_data%lagsatex)
        call sqlite_insert_row_add_real(wb_data%wet_evap)
        call sqlite_insert_row_add_real(wb_data%wet_out)
        call sqlite_insert_row_add_real(wb_data%wet_stor)
        call sqlite_insert_row_add_text(plant_cov)
        call sqlite_insert_row_add_text(mgt_ops)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_hru_wb
      
      !> Insert HRU nutrient balance data into SQLite
      subroutine sqlite_insert_hru_nb(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, nb_data, plant_cov, mgt_ops)
        use output_landscape_module, only : output_nutbal
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name, plant_cov, mgt_ops
        type(output_nutbal), intent(in) :: nb_data
        character(len=32) :: table_name
        character(len=16) :: col_names(29)
        
        if (.not. sqlite_initialized) return
        
        table_name = "hru_nb_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "grazn", "grazp", "lab_min_p", "act_sta_p", "fertn", "fertp", "fixn", "denit", &
                     "act_nit_n", "act_sta_n", "org_lab_p", "rsd_nitorg_n", "rsd_laborg_p", "no3atmo", &
                     "nh4atmo", "nuptake", "puptake", "gwsoiln", "gwsoilp", "plant_cov", "mgt_ops", "null1"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(nb_data%grazn)
        call sqlite_insert_row_add_real(nb_data%grazp)
        call sqlite_insert_row_add_real(nb_data%lab_min_p)
        call sqlite_insert_row_add_real(nb_data%act_sta_p)
        call sqlite_insert_row_add_real(nb_data%fertn)
        call sqlite_insert_row_add_real(nb_data%fertp)
        call sqlite_insert_row_add_real(nb_data%fixn)
        call sqlite_insert_row_add_real(nb_data%denit)
        call sqlite_insert_row_add_real(nb_data%act_nit_n)
        call sqlite_insert_row_add_real(nb_data%act_sta_n)
        call sqlite_insert_row_add_real(nb_data%org_lab_p)
        call sqlite_insert_row_add_real(nb_data%rsd_nitorg_n)
        call sqlite_insert_row_add_real(nb_data%rsd_laborg_p)
        call sqlite_insert_row_add_real(nb_data%no3atmo)
        call sqlite_insert_row_add_real(nb_data%nh4atmo)
        call sqlite_insert_row_add_real(nb_data%nuptake)
        call sqlite_insert_row_add_real(nb_data%puptake)
        call sqlite_insert_row_add_real(nb_data%gwsoiln)
        call sqlite_insert_row_add_real(nb_data%gwsoilp)
        call sqlite_insert_row_add_text(plant_cov)
        call sqlite_insert_row_add_text(mgt_ops)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_hru_nb
      
      !> Insert HRU losses data into SQLite
      subroutine sqlite_insert_hru_ls(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, ls_data, plant_cov, mgt_ops)
        use output_landscape_module, only : output_losses
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name, plant_cov, mgt_ops
        type(output_losses), intent(in) :: ls_data
        character(len=32) :: table_name
        character(len=16) :: col_names(21)
        
        if (.not. sqlite_initialized) return
        
        table_name = "hru_ls_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "sedyld", "sedorgn", "sedorgp", "surqno3", "latno3", "surqsolp", "usle", &
                     "sedminp", "tileno3", "lchlabp", "tilelabp", "satexn", "plant_cov", "mgt_ops"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(ls_data%sedyld)
        call sqlite_insert_row_add_real(ls_data%sedorgn)
        call sqlite_insert_row_add_real(ls_data%sedorgp)
        call sqlite_insert_row_add_real(ls_data%surqno3)
        call sqlite_insert_row_add_real(ls_data%latno3)
        call sqlite_insert_row_add_real(ls_data%surqsolp)
        call sqlite_insert_row_add_real(ls_data%usle)
        call sqlite_insert_row_add_real(ls_data%sedminp)
        call sqlite_insert_row_add_real(ls_data%tileno3)
        call sqlite_insert_row_add_real(ls_data%lchlabp)
        call sqlite_insert_row_add_real(ls_data%tilelabp)
        call sqlite_insert_row_add_real(ls_data%satexn)
        call sqlite_insert_row_add_text(plant_cov)
        call sqlite_insert_row_add_text(mgt_ops)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_hru_ls
      
      !> Insert HRU plant weather data into SQLite
      subroutine sqlite_insert_hru_pw(period, jday, mon, day_mo, yrc, unit_id, gis_id, name, pw_data, plant_cov, mgt_ops)
        use output_landscape_module, only : output_plantweather
        
        character(len=*), intent(in) :: period
        integer, intent(in) :: jday, mon, day_mo, yrc, unit_id
        integer*8, intent(in) :: gis_id
        character(len=*), intent(in) :: name, plant_cov, mgt_ops
        type(output_plantweather), intent(in) :: pw_data
        character(len=32) :: table_name
        character(len=16) :: col_names(34)
        
        if (.not. sqlite_initialized) return
        
        table_name = "hru_pw_" // trim(period)
        
        col_names = [character(len=16) :: "jday", "mon", "day", "yr", "unit", "gis_id", "name", &
                     "lai", "bioms", "yield", "residue", "sol_tmp", "strsw", "strsa", "strstmp", &
                     "strsn", "strsp", "strss", "nplnt", "percn", "pplnt", "tmx", "tmn", "tmpav", &
                     "solrad", "wndspd", "rhum", "phubase0", "lai_max", "bm_max", "bm_grow", "c_gro", &
                     "plant_cov", "mgt_ops"]
        
        call sqlite_insert_row_start_named(trim(table_name), col_names)
        call sqlite_insert_row_add_int(jday)
        call sqlite_insert_row_add_int(mon)
        call sqlite_insert_row_add_int(day_mo)
        call sqlite_insert_row_add_int(yrc)
        call sqlite_insert_row_add_int(unit_id)
        call sqlite_insert_row_add_int8(gis_id)
        call sqlite_insert_row_add_text(name)
        call sqlite_insert_row_add_real(pw_data%lai)
        call sqlite_insert_row_add_real(pw_data%bioms)
        call sqlite_insert_row_add_real(pw_data%yield)
        call sqlite_insert_row_add_real(pw_data%residue)
        call sqlite_insert_row_add_real(pw_data%sol_tmp)
        call sqlite_insert_row_add_real(pw_data%strsw)
        call sqlite_insert_row_add_real(pw_data%strsa)
        call sqlite_insert_row_add_real(pw_data%strstmp)
        call sqlite_insert_row_add_real(pw_data%strsn)
        call sqlite_insert_row_add_real(pw_data%strsp)
        call sqlite_insert_row_add_real(pw_data%strss)
        call sqlite_insert_row_add_real(pw_data%nplnt)
        call sqlite_insert_row_add_real(pw_data%percn)
        call sqlite_insert_row_add_real(pw_data%pplnt)
        call sqlite_insert_row_add_real(pw_data%tmx)
        call sqlite_insert_row_add_real(pw_data%tmn)
        call sqlite_insert_row_add_real(pw_data%tmpav)
        call sqlite_insert_row_add_real(pw_data%solrad)
        call sqlite_insert_row_add_real(pw_data%wndspd)
        call sqlite_insert_row_add_real(pw_data%rhum)
        call sqlite_insert_row_add_real(pw_data%phubase0)
        call sqlite_insert_row_add_real(pw_data%lai_max)
        call sqlite_insert_row_add_real(pw_data%bm_max)
        call sqlite_insert_row_add_real(pw_data%bm_grow)
        call sqlite_insert_row_add_real(pw_data%c_gro)
        call sqlite_insert_row_add_text(plant_cov)
        call sqlite_insert_row_add_text(mgt_ops)
        call sqlite_insert_row_execute()
        
      end subroutine sqlite_insert_hru_pw
      
      end module sqlite_output_module
      
      !> Finalize SQLite output - called at end of simulation
      !! This is a standalone subroutine to be called from main program
      subroutine sqlite_finalize_output()
        use sqlite_output_module, only : sqlite_finalize
        
        implicit none
        
        call sqlite_finalize()
        write(*,*) "SQLite output database finalized."
        
      end subroutine sqlite_finalize_output
#endif
