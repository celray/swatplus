#ifdef SQLITE_ENABLED
      module sqlite3_module
      !! SQLite3 C API bindings using ISO_C_BINDING
      !! Works with gfortran, ifx, and ifort
      !! Requires linking with -lsqlite3 (Unix) or sqlite3.lib (Windows)
      
      use, intrinsic :: iso_c_binding
      implicit none

      private

      ! Public constants
      public :: SQLITE_OK, SQLITE_ERROR, SQLITE_ROW, SQLITE_DONE
      public :: SQLITE_INTEGER, SQLITE_FLOAT, SQLITE_TEXT, SQLITE_BLOB, SQLITE_NULL

      ! Public types
      public :: sqlite3_db, sqlite3_stmt_ptr

      ! Public procedures
      public :: sqlite3_open, sqlite3_close
      public :: sqlite3_exec, sqlite3_errmsg_ptr, sqlite3_errmsg
      public :: sqlite3_prepare_v2, sqlite3_step, sqlite3_finalize, sqlite3_reset, sqlite3_clear_bindings
      public :: sqlite3_bind_int, sqlite3_bind_int64, sqlite3_bind_double, sqlite3_bind_text, sqlite3_bind_null
      public :: sqlite3_column_int, sqlite3_column_double, sqlite3_column_text_ptr
      public :: sqlite3_last_insert_rowid
      public :: fortran_to_c_string, c_to_fortran_string

      ! SQLite result codes
      integer(c_int), parameter :: SQLITE_OK = 0
      integer(c_int), parameter :: SQLITE_ERROR = 1
      integer(c_int), parameter :: SQLITE_ROW = 100
      integer(c_int), parameter :: SQLITE_DONE = 101

      ! SQLite data types
      integer(c_int), parameter :: SQLITE_INTEGER = 1
      integer(c_int), parameter :: SQLITE_FLOAT = 2
      integer(c_int), parameter :: SQLITE_TEXT = 3
      integer(c_int), parameter :: SQLITE_BLOB = 4
      integer(c_int), parameter :: SQLITE_NULL = 5

      ! Type aliases for clarity
      type :: sqlite3_db
        type(c_ptr) :: ptr = c_null_ptr
      end type sqlite3_db

      type :: sqlite3_stmt_ptr
        type(c_ptr) :: ptr = c_null_ptr
      end type sqlite3_stmt_ptr

      ! C function interfaces
      interface

        !! Open a database connection
        function c_sqlite3_open(filename, ppDb) bind(C, name='sqlite3_open')
          import :: c_ptr, c_char, c_int
          character(kind=c_char), intent(in) :: filename(*)
          type(c_ptr), intent(out) :: ppDb
          integer(c_int) :: c_sqlite3_open
        end function c_sqlite3_open

        !! Close a database connection
        function c_sqlite3_close(db) bind(C, name='sqlite3_close')
          import :: c_ptr, c_int
          type(c_ptr), value :: db
          integer(c_int) :: c_sqlite3_close
        end function c_sqlite3_close

        !! Execute SQL directly (for DDL and simple statements)
        function c_sqlite3_exec(db, sql, callback, arg, errmsg) bind(C, name='sqlite3_exec')
          import :: c_ptr, c_char, c_int, c_funptr
          type(c_ptr), value :: db
          character(kind=c_char), intent(in) :: sql(*)
          type(c_funptr), value :: callback
          type(c_ptr), value :: arg
          type(c_ptr) :: errmsg
          integer(c_int) :: c_sqlite3_exec
        end function c_sqlite3_exec

        !! Get error message pointer
        function c_sqlite3_errmsg(db) bind(C, name='sqlite3_errmsg')
          import :: c_ptr
          type(c_ptr), value :: db
          type(c_ptr) :: c_sqlite3_errmsg
        end function c_sqlite3_errmsg

        !! Prepare a statement
        function c_sqlite3_prepare_v2(db, sql, nByte, ppStmt, pzTail) bind(C, name='sqlite3_prepare_v2')
          import :: c_ptr, c_char, c_int
          type(c_ptr), value :: db
          character(kind=c_char), intent(in) :: sql(*)
          integer(c_int), value :: nByte
          type(c_ptr), intent(out) :: ppStmt
          type(c_ptr), intent(out) :: pzTail
          integer(c_int) :: c_sqlite3_prepare_v2
        end function c_sqlite3_prepare_v2

        !! Step through a prepared statement
        function c_sqlite3_step(pStmt) bind(C, name='sqlite3_step')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int) :: c_sqlite3_step
        end function c_sqlite3_step

        !! Finalize (destroy) a prepared statement
        function c_sqlite3_finalize(pStmt) bind(C, name='sqlite3_finalize')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int) :: c_sqlite3_finalize
        end function c_sqlite3_finalize

        !! Reset a prepared statement for re-use
        function c_sqlite3_reset(pStmt) bind(C, name='sqlite3_reset')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int) :: c_sqlite3_reset
        end function c_sqlite3_reset

        !! Clear all bindings on a prepared statement
        function c_sqlite3_clear_bindings(pStmt) bind(C, name='sqlite3_clear_bindings')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int) :: c_sqlite3_clear_bindings
        end function c_sqlite3_clear_bindings

        !! Bind integer value
        function c_sqlite3_bind_int(pStmt, idx, val) bind(C, name='sqlite3_bind_int')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int), value :: idx
          integer(c_int), value :: val
          integer(c_int) :: c_sqlite3_bind_int
        end function c_sqlite3_bind_int

        !! Bind 64-bit integer value
        function c_sqlite3_bind_int64(pStmt, idx, val) bind(C, name='sqlite3_bind_int64')
          import :: c_ptr, c_int, c_int64_t
          type(c_ptr), value :: pStmt
          integer(c_int), value :: idx
          integer(c_int64_t), value :: val
          integer(c_int) :: c_sqlite3_bind_int64
        end function c_sqlite3_bind_int64

        !! Bind double value
        function c_sqlite3_bind_double(pStmt, idx, val) bind(C, name='sqlite3_bind_double')
          import :: c_ptr, c_int, c_double
          type(c_ptr), value :: pStmt
          integer(c_int), value :: idx
          real(c_double), value :: val
          integer(c_int) :: c_sqlite3_bind_double
        end function c_sqlite3_bind_double

        !! Bind text value
        function c_sqlite3_bind_text(pStmt, idx, val, nByte, destructor) bind(C, name='sqlite3_bind_text')
          import :: c_ptr, c_char, c_int, c_funptr
          type(c_ptr), value :: pStmt
          integer(c_int), value :: idx
          character(kind=c_char), intent(in) :: val(*)
          integer(c_int), value :: nByte
          type(c_funptr), value :: destructor
          integer(c_int) :: c_sqlite3_bind_text
        end function c_sqlite3_bind_text

        !! Bind NULL value
        function c_sqlite3_bind_null(pStmt, idx) bind(C, name='sqlite3_bind_null')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int), value :: idx
          integer(c_int) :: c_sqlite3_bind_null
        end function c_sqlite3_bind_null

        !! Get integer column value
        function c_sqlite3_column_int(pStmt, iCol) bind(C, name='sqlite3_column_int')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int), value :: iCol
          integer(c_int) :: c_sqlite3_column_int
        end function c_sqlite3_column_int

        !! Get double column value
        function c_sqlite3_column_double(pStmt, iCol) bind(C, name='sqlite3_column_double')
          import :: c_ptr, c_int, c_double
          type(c_ptr), value :: pStmt
          integer(c_int), value :: iCol
          real(c_double) :: c_sqlite3_column_double
        end function c_sqlite3_column_double

        !! Get text column value (returns pointer)
        function c_sqlite3_column_text(pStmt, iCol) bind(C, name='sqlite3_column_text')
          import :: c_ptr, c_int
          type(c_ptr), value :: pStmt
          integer(c_int), value :: iCol
          type(c_ptr) :: c_sqlite3_column_text
        end function c_sqlite3_column_text

        !! Get last insert rowid
        function c_sqlite3_last_insert_rowid(db) bind(C, name='sqlite3_last_insert_rowid')
          import :: c_ptr, c_int64_t
          type(c_ptr), value :: db
          integer(c_int64_t) :: c_sqlite3_last_insert_rowid
        end function c_sqlite3_last_insert_rowid

        !! C strlen for string conversion
        function c_strlen(s) bind(C, name='strlen')
          import :: c_ptr, c_size_t
          type(c_ptr), value :: s
          integer(c_size_t) :: c_strlen
        end function c_strlen

      end interface

      ! SQLITE_TRANSIENT constant (-1 cast to destructor function pointer)
      ! This tells SQLite to make its own copy of the string
      type(c_funptr), parameter :: SQLITE_TRANSIENT = transfer(-1_c_intptr_t, c_null_funptr)

      contains

      !> Convert Fortran string to C string (null-terminated) - returns allocatable array
      function fortran_to_c_string(fstring) result(cstring)
        character(len=*), intent(in) :: fstring
        character(kind=c_char), allocatable :: cstring(:)
        integer :: i, n

        n = len_trim(fstring)
        allocate(cstring(n + 1))
        do i = 1, n
          cstring(i) = fstring(i:i)
        end do
        cstring(n + 1) = c_null_char
      end function fortran_to_c_string

      !> Convert C string pointer to Fortran string
      function c_to_fortran_string(cptr) result(fstring)
        type(c_ptr), intent(in) :: cptr
        character(len=256) :: fstring
        character(kind=c_char), pointer :: carr(:)
        integer :: i, n

        fstring = ''
        if (.not. c_associated(cptr)) return

        n = int(c_strlen(cptr))
        if (n > 256) n = 256
        if (n == 0) return

        call c_f_pointer(cptr, carr, [n])
        do i = 1, n
          fstring(i:i) = carr(i)
        end do
      end function c_to_fortran_string

      !> Open a SQLite database
      function sqlite3_open(filename, db) result(rc)
        character(len=*), intent(in) :: filename
        type(sqlite3_db), intent(out) :: db
        integer(c_int) :: rc
        character(kind=c_char), allocatable :: c_filename(:)

        c_filename = fortran_to_c_string(filename)
        rc = c_sqlite3_open(c_filename, db%ptr)
      end function sqlite3_open

      !> Close a SQLite database
      function sqlite3_close(db) result(rc)
        type(sqlite3_db), intent(inout) :: db
        integer(c_int) :: rc

        rc = c_sqlite3_close(db%ptr)
        db%ptr = c_null_ptr
      end function sqlite3_close

      !> Execute SQL directly (no callback, no results)
      function sqlite3_exec(db, sql) result(rc)
        type(sqlite3_db), intent(in) :: db
        character(len=*), intent(in) :: sql
        integer(c_int) :: rc
        type(c_ptr) :: errmsg
        character(kind=c_char), allocatable :: c_sql(:)

        c_sql = fortran_to_c_string(sql)
        rc = c_sqlite3_exec(db%ptr, c_sql, c_null_funptr, c_null_ptr, errmsg)
      end function sqlite3_exec

      !> Get error message pointer (for advanced use)
      function sqlite3_errmsg_ptr(db) result(ptr)
        type(sqlite3_db), intent(in) :: db
        type(c_ptr) :: ptr

        ptr = c_sqlite3_errmsg(db%ptr)
      end function sqlite3_errmsg_ptr

      !> Get error message as Fortran string
      function sqlite3_errmsg(db) result(msg)
        type(sqlite3_db), intent(in) :: db
        character(len=256) :: msg
        type(c_ptr) :: ptr

        ptr = c_sqlite3_errmsg(db%ptr)
        msg = c_to_fortran_string(ptr)
      end function sqlite3_errmsg

      !> Prepare a SQL statement
      function sqlite3_prepare_v2(db, sql, stmt) result(rc)
        type(sqlite3_db), intent(in) :: db
        character(len=*), intent(in) :: sql
        type(sqlite3_stmt_ptr), intent(out) :: stmt
        integer(c_int) :: rc
        type(c_ptr) :: tail
        character(kind=c_char), allocatable :: c_sql(:)

        c_sql = fortran_to_c_string(sql)
        rc = c_sqlite3_prepare_v2(db%ptr, c_sql, -1_c_int, stmt%ptr, tail)
      end function sqlite3_prepare_v2

      !> Step through a prepared statement
      function sqlite3_step(stmt) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer(c_int) :: rc

        rc = c_sqlite3_step(stmt%ptr)
      end function sqlite3_step

      !> Finalize a prepared statement
      function sqlite3_finalize(stmt) result(rc)
        type(sqlite3_stmt_ptr), intent(inout) :: stmt
        integer(c_int) :: rc

        rc = c_sqlite3_finalize(stmt%ptr)
        stmt%ptr = c_null_ptr
      end function sqlite3_finalize

      !> Reset a prepared statement
      function sqlite3_reset(stmt) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer(c_int) :: rc
        
        rc = c_sqlite3_reset(stmt%ptr)
      end function sqlite3_reset

      !> Clear all bindings on a prepared statement
      function sqlite3_clear_bindings(stmt) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer(c_int) :: rc
        
        rc = c_sqlite3_clear_bindings(stmt%ptr)
      end function sqlite3_clear_bindings

      !> Bind integer value (1-indexed)
      function sqlite3_bind_int(stmt, idx, val) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: idx
        integer, intent(in) :: val
        integer(c_int) :: rc

        rc = c_sqlite3_bind_int(stmt%ptr, int(idx, c_int), int(val, c_int))
      end function sqlite3_bind_int

      !> Bind 64-bit integer value (1-indexed)
      function sqlite3_bind_int64(stmt, idx, val) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: idx
        integer*8, intent(in) :: val
        integer(c_int) :: rc

        rc = c_sqlite3_bind_int64(stmt%ptr, int(idx, c_int), int(val, c_int64_t))
      end function sqlite3_bind_int64

      !> Bind double/real value (1-indexed)
      function sqlite3_bind_double(stmt, idx, val) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: idx
        real(8), intent(in) :: val
        integer(c_int) :: rc

        rc = c_sqlite3_bind_double(stmt%ptr, int(idx, c_int), real(val, c_double))
      end function sqlite3_bind_double

      !> Bind text value (1-indexed)
      function sqlite3_bind_text(stmt, idx, val) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: idx
        character(len=*), intent(in) :: val
        integer(c_int) :: rc
        character(kind=c_char), allocatable :: cstr(:)

        cstr = fortran_to_c_string(val)
        rc = c_sqlite3_bind_text(stmt%ptr, int(idx, c_int), cstr, int(len_trim(val), c_int), SQLITE_TRANSIENT)
      end function sqlite3_bind_text

      !> Bind NULL value (1-indexed)
      function sqlite3_bind_null(stmt, idx) result(rc)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: idx
        integer(c_int) :: rc

        rc = c_sqlite3_bind_null(stmt%ptr, int(idx, c_int))
      end function sqlite3_bind_null

      !> Get integer column value (0-indexed)
      function sqlite3_column_int(stmt, iCol) result(val)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: iCol
        integer :: val

        val = int(c_sqlite3_column_int(stmt%ptr, int(iCol, c_int)))
      end function sqlite3_column_int

      !> Get double column value (0-indexed)
      function sqlite3_column_double(stmt, iCol) result(val)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: iCol
        real(8) :: val

        val = real(c_sqlite3_column_double(stmt%ptr, int(iCol, c_int)), 8)
      end function sqlite3_column_double

      !> Get text column pointer (0-indexed) - use c_to_fortran_string to convert
      function sqlite3_column_text_ptr(stmt, iCol) result(ptr)
        type(sqlite3_stmt_ptr), intent(in) :: stmt
        integer, intent(in) :: iCol
        type(c_ptr) :: ptr

        ptr = c_sqlite3_column_text(stmt%ptr, int(iCol, c_int))
      end function sqlite3_column_text_ptr

      !> Get last insert rowid
      function sqlite3_last_insert_rowid(db) result(rowid)
        type(sqlite3_db), intent(in) :: db
        integer(8) :: rowid

        rowid = int(c_sqlite3_last_insert_rowid(db%ptr), 8)
      end function sqlite3_last_insert_rowid

      end module sqlite3_module
#endif
