subroutine cli_ncdf_meas
    
    ! This subroutine reads NetCDF climate data for SWAT+ simulation

    use climate_module
    use maximum_data_module
    use time_module
    use input_file_module
    
    implicit none
    
    ! Use C-style NetCDF interface (works with Intel compilers)
    include 'netcdf.inc'
    
    ! NetCDF variables
    integer :: ncid, varid, dimid
    integer :: status
    character(len=256) :: ncdf_file
    character(len=256) :: version_string
    
    ! NetCDF dimensions
    integer :: time_dimid, lat_dimid, lon_dimid
    integer :: ntime, nlat, nlon
    
    ! NetCDF variable IDs
    integer :: time_varid, lat_varid, lon_varid
    integer :: pcp_varid, tmin_varid, tmax_varid, slr_varid, hmd_varid, wnd_varid
    
    ! Arrays for reading ALL data
    real, dimension(:), allocatable :: time_vals, lat_vals, lon_vals
    real, dimension(:,:,:), allocatable :: pcp_data, tmin_data, tmax_data
    real, dimension(:,:,:), allocatable :: slr_data, hmd_data, wnd_data
    
    ! Variables for finding closest grid point
    integer :: target_lat_idx, target_lon_idx
    real :: target_lat, target_lon, min_dist, dist
    integer :: ilat, ilon
    
    ! Variables for date calculation
    integer :: days_since_ref, year, month, day, days_in_months(12)
    integer :: ref_year, ref_month, ref_day
    character(len=256) :: time_units
    
    ! Loop counters and diagnostics
    integer :: istat, itime, iyear, iday, i, iwst
    integer :: start_year, start_day, current_year, days_in_year
    integer :: actual_year  ! For leap year calculation
    integer :: days_in_year_loop ! For proper leap year counting
    integer :: ndims, dimids(NF_MAX_VAR_DIMS), natts, xtype
    character(len=NF_MAX_NAME) :: var_name, dim_name
    integer :: dim_len
    
    ! Variables for monthly statistics calculation
    integer :: mo, day_mo
    logical :: exists
    
    version_string = nf_inq_libvers()
    if (version_string(len_trim(version_string):len_trim(version_string)) == '$') then
        version_string(len_trim(version_string):len_trim(version_string)) = ' '
    endif
    write (*,*) "reading data using netcdf version ", trim(version_string)
    write (9003,*) "reading data using netcdf version ", trim(version_string)

    ! Get NetCDF filename based on precipitation path from file.cio
    if (in_path_pcp%pcp == "null" .or. trim(in_path_pcp%pcp) == " ") then
        ncdf_file = ""
        write(*,*) "! error: No NetCDF file specified in pcp path in 'file.cio'"
        write (9003,*) "! error: No NetCDF file specified in pcp path in 'file.cio'"
        stop
    else
        ! check if the file exists
        ncdf_file = TRIM(ADJUSTL(in_path_pcp%pcp))
        inquire(file=trim(ncdf_file), exist=exists)
        if (.not. exists) then
            write(*,*) "! error: NetCDF file does not exist at", trim(ncdf_file)
            write (9003,*) "! error: NetCDF file does not exist at", trim(ncdf_file)
            stop
        end if
    endif
    
    ! Open NetCDF file
    status = nf_open(trim(ncdf_file), NF_NOWRITE, ncid)
    if (status /= NF_NOERR) then
        write (*,*) "ERROR: Cannot open NetCDF file: ", trim(ncdf_file)
        write (*,*) "NetCDF Error: ", nf_strerror(status)
        write (9003,*) "ERROR: Cannot open NetCDF file: ", trim(ncdf_file)
        stop
    endif
    
    ! Get dimensions for gridded data
    status = nf_inq_dimid(ncid, "time", time_dimid)
    if (status /= NF_NOERR) then
        write (*,*) "ERROR: Cannot find 'time' dimension"
        write (9003,*) "ERROR: Cannot find 'time' dimension"
        stop
    endif
    
    status = nf_inq_dimid(ncid, "lat", lat_dimid)
    if (status /= NF_NOERR) then
        write (*,*) "ERROR: Cannot find 'lat' dimension"
        write (9003,*) "ERROR: Cannot find 'lat' dimension"
        stop
    endif
    
    status = nf_inq_dimid(ncid, "lon", lon_dimid)
    if (status /= NF_NOERR) then
        write (*,*) "ERROR: Cannot find 'lon' dimension"
        write (9003,*) "ERROR: Cannot find 'lon' dimension"
        stop
    endif
    
    ! Get dimension sizes
    status = nf_inq_dimlen(ncid, time_dimid, ntime)
    status = nf_inq_dimlen(ncid, lat_dimid, nlat)
    status = nf_inq_dimlen(ncid, lon_dimid, nlon)
    
    ! Check variable dimension order before allocating arrays
    status = nf_inq_varid(ncid, "pcp", pcp_varid)
    if (status == NF_NOERR) then
        status = nf_inq_var(ncid, pcp_varid, var_name, xtype, ndims, dimids, natts)
        do i = 1, ndims
            status = nf_inq_dim(ncid, dimids(i), dim_name, dim_len)
        end do
    endif
    
    ! Allocate arrays - actual storage order is (lon, lat, time)
    allocate(time_vals(ntime), lat_vals(nlat), lon_vals(nlon))
    allocate(pcp_data(nlon, nlat, ntime), tmin_data(nlon, nlat, ntime), tmax_data(nlon, nlat, ntime))
    allocate(slr_data(nlon, nlat, ntime), hmd_data(nlon, nlat, ntime), wnd_data(nlon, nlat, ntime))
    
    ! Read coordinate data
    status = nf_inq_varid(ncid, "lat", lat_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, lat_varid, lat_vals)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading lat data: ", nf_strerror(status)
            write (9003,*) "Error reading lat data: ", nf_strerror(status)
            stop
        endif
    endif
    
    status = nf_inq_varid(ncid, "lon", lon_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, lon_varid, lon_vals)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading lon data: ", nf_strerror(status)
            write (9003,*) "Error reading lon data: ", nf_strerror(status)
            stop
        endif
    endif
    
    ! Skip to 100 continue
    goto 100
    
100 continue
    
    ! Read time data to get the date for first time step
    status = nf_inq_varid(ncid, "time", time_varid)
    if (status == NF_NOERR) then
        ! Read time units attribute to get reference date
        status = nf_get_att_text(ncid, time_varid, "units", time_units)
        if (status /= NF_NOERR) then
            write (*,*) "! warning: Cannot read time units attribute, assuming days since 1970-01-01"
            write (9003,*) "! warning: Cannot read time units attribute, assuming days since 1970-01-01"
            time_units = "days since 1970-01-01 00:00:00"
        endif
        
        ! Parse reference date from time units string
        ! Expected format: "days since YYYY-MM-DD HH:MM:SS" or "days since YYYY-MM-DD"
        call parse_time_units(time_units, ref_year, ref_month, ref_day)
        
        status = nf_get_var_real(ncid, time_varid, time_vals)
        if (status == NF_NOERR) then
            if (ntime >= 1) then
                ! Time is in "days since ref_date"
                days_since_ref = int(time_vals(1))
                
                ! Calculate actual date from reference date + days offset
                call add_days_to_date(ref_year, ref_month, ref_day, days_since_ref, year, month, day)
            endif
        else
            write (*,*) "Error reading time data: ", nf_strerror(status)
            write (9003,*) "Error reading time data: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No time variable found in NetCDF"
        write (9003,*) "WARNING: No time variable found in NetCDF"
        stop
    endif
    
    ! Precipitation
    status = nf_inq_varid(ncid, "pcp", pcp_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, pcp_varid, pcp_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading precipitation: ", nf_strerror(status)
            write (9003,*) "Error reading precipitation: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No pcp variable found in NetCDF"
        write (9003,*) "WARNING: No pcp variable found in NetCDF"
    endif
    
    ! Temperature minimum
    status = nf_inq_varid(ncid, "tmin", tmin_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, tmin_varid, tmin_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading tmin: ", nf_strerror(status)
            write (9003,*) "Error reading tmin: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No tmin variable found in NetCDF"
        write (9003,*) "WARNING: No tmin variable found in NetCDF"
    endif
    
    ! Temperature maximum  
    status = nf_inq_varid(ncid, "tmax", tmax_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, tmax_varid, tmax_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading tmax: ", nf_strerror(status)
            write (9003,*) "Error reading tmax: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No tmax variable found in NetCDF"
        write (9003,*) "WARNING: No tmax variable found in NetCDF"
    endif
    
    ! Solar radiation
    status = nf_inq_varid(ncid, "slr", slr_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, slr_varid, slr_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading slr: ", nf_strerror(status)
            write (9003,*) "Error reading slr: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No slr variable found in NetCDF"
        write (9003,*) "WARNING: No slr variable found in NetCDF"
    endif
    
    ! Relative humidity
    status = nf_inq_varid(ncid, "hmd", hmd_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, hmd_varid, hmd_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading hmd: ", nf_strerror(status)
            write (9003,*) "Error reading hmd: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No hmd variable found in NetCDF"
        write (9003,*) "WARNING: No hmd variable found in NetCDF"
    endif
    
    ! Wind speed
    status = nf_inq_varid(ncid, "wnd", wnd_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, wnd_varid, wnd_data)
        if (status /= NF_NOERR) then
            write (*,*) "Error reading wnd: ", nf_strerror(status)
            write (9003,*) "Error reading wnd: ", nf_strerror(status)
            stop
        endif
    else
        write (*,*) "WARNING: No wnd variable found in NetCDF"
        write (9003,*) "WARNING: No wnd variable found in NetCDF"
    endif

    ! allocate and populate climate arrays
    allocate (pcp(0:db_mx%wst))
    allocate (pcp_n(db_mx%wst))
    allocate (tmp(0:db_mx%wst))
    allocate (tmp_n(db_mx%wst))
    allocate (slr(0:db_mx%wst))
    allocate (slr_n(db_mx%wst))
    allocate (hmd(0:db_mx%wst))
    allocate (hmd_n(db_mx%wst))
    allocate (wnd(0:db_mx%wst))
    allocate (wnd_n(db_mx%wst))
    db_mx%pcpfiles = db_mx%wst
    db_mx%tmpfiles = db_mx%wst  
    db_mx%slrfiles = db_mx%wst
    db_mx%rhfiles = db_mx%wst
    db_mx%wndfiles = db_mx%wst
    
    ! Set the climate file indices for each station based on station index
    ! This is needed because NetCDF data is indexed by station, unlike traditional
    ! files which use the search() function to map file names to indices
    do i = 1, db_mx%wst
        wst(i)%wco%pgage = i
        wst(i)%wco%tgage = i
        wst(i)%wco%sgage = i
        wst(i)%wco%hgage = i
        wst(i)%wco%wgage = i
    end do
    
    do iwst = 1, db_mx%wst
        
        ! Find closest grid point to this station
        min_dist = 999999.0
        target_lat_idx = 1
        target_lon_idx = 1
        
        do ilat = 1, nlat
            do ilon = 1, nlon
                dist = sqrt((lat_vals(ilat) - wst(iwst)%lat)**2 + (lon_vals(ilon) - wst(iwst)%lon)**2)
                if (dist < min_dist) then
                    min_dist = dist
                    target_lat_idx = ilat
                    target_lon_idx = ilon
                endif
            end do
        end do
        
        ! Set station names as "filenames" 
        pcp_n(iwst) = wst(iwst)%name
        tmp_n(iwst) = wst(iwst)%name
        slr_n(iwst) = wst(iwst)%name 
        hmd_n(iwst) = wst(iwst)%name
        wnd_n(iwst) = wst(iwst)%name
        
        ! Populate station metadata from wst and NetCDF grid
        ! CRITICAL: Use station coordinates (not NetCDF grid) to match traditional method
        pcp(iwst)%filename = wst(iwst)%name
        pcp(iwst)%lat = wst(iwst)%lat      ! Use station coordinates, not NetCDF grid
        pcp(iwst)%long = wst(iwst)%lon     ! Use station coordinates, not NetCDF grid
        pcp(iwst)%elev = wst(iwst)%elev    ! Use station elevation
        
        ! Copy metadata to other climate arrays
        tmp(iwst)%filename = pcp(iwst)%filename
        tmp(iwst)%lat = pcp(iwst)%lat
        tmp(iwst)%long = pcp(iwst)%long
        tmp(iwst)%elev = pcp(iwst)%elev
        
        slr(iwst)%filename = pcp(iwst)%filename
        slr(iwst)%lat = pcp(iwst)%lat
        slr(iwst)%long = pcp(iwst)%long
        slr(iwst)%elev = pcp(iwst)%elev
        
        hmd(iwst)%filename = pcp(iwst)%filename
        hmd(iwst)%lat = pcp(iwst)%lat
        hmd(iwst)%long = pcp(iwst)%long
        hmd(iwst)%elev = pcp(iwst)%elev
        
        wnd(iwst)%filename = pcp(iwst)%filename
        wnd(iwst)%lat = pcp(iwst)%lat
        wnd(iwst)%long = pcp(iwst)%long
        wnd(iwst)%elev = pcp(iwst)%elev
        
    end do

    ! Step 2: Populate the actual time series data for each station
    ! Similar to how cli_pmeas.f90 populates pcp(i)%ts(day, year)
    
    do iwst = 1, db_mx%wst
        
        ! Find closest grid point to this station (recompute for each station)
        min_dist = 999999.0
        target_lat_idx = 1
        target_lon_idx = 1
        
        do ilat = 1, nlat
            do ilon = 1, nlon
                dist = sqrt((lat_vals(ilat) - wst(iwst)%lat)**2 + (lon_vals(ilon) - wst(iwst)%lon)**2)
                if (dist < min_dist) then
                    min_dist = dist
                    target_lat_idx = ilat
                    target_lon_idx = ilon
                endif
            end do
        end do
        
        ! Calculate number of years from NetCDF time dimension
        ! Assume daily data, so nyears = ntime / 365 (approximately)
        ! For more accurate calculation, we'd need to parse the actual dates
        pcp(iwst)%nbyr = max(1, ntime / 365)
        tmp(iwst)%nbyr = pcp(iwst)%nbyr
        slr(iwst)%nbyr = pcp(iwst)%nbyr
        hmd(iwst)%nbyr = pcp(iwst)%nbyr
        wnd(iwst)%nbyr = pcp(iwst)%nbyr
        
        ! Set timestep (0 = daily, >0 = sub-daily)
        pcp(iwst)%tstep = 0  ! Assume daily for now
        tmp(iwst)%tstep = 0
        slr(iwst)%tstep = 0
        hmd(iwst)%tstep = 0
        wnd(iwst)%tstep = 0
        
        ! Initialize days_gen counter (like traditional method)
        pcp(iwst)%days_gen = 0
        tmp(iwst)%days_gen = 0
        slr(iwst)%days_gen = 0
        hmd(iwst)%days_gen = 0
        wnd(iwst)%days_gen = 0
        
        ! Set start and end years (approximate)
        pcp(iwst)%start_yr = time%yrc
        pcp(iwst)%end_yr = time%yrc + pcp(iwst)%nbyr - 1
        pcp(iwst)%start_day = 1
        
        ! Calculate proper end_day based on leap year for last year (like traditional method)
        actual_year = pcp(iwst)%end_yr
        ! Use proper leap year calculation (same as time_control.f90)
        if (Mod(actual_year,4) == 0) then
          if (Mod(actual_year,100) == 0) then
            if (Mod(actual_year,400) == 0) then
              pcp(iwst)%end_day = 366  ! Last day of leap year
            else
              pcp(iwst)%end_day = 365  ! Century year, not leap
            end if
          else
            pcp(iwst)%end_day = 366  ! Last day of leap year
          end if
        else
          pcp(iwst)%end_day = 365  ! Last day of non-leap year
        endif
        
        ! Calculate yrs_start correctly (like traditional method)
        if (pcp(iwst)%start_yr > time%yrc) then
            pcp(iwst)%yrs_start = pcp(iwst)%start_yr - time%yrc
        else
            pcp(iwst)%yrs_start = 0
        end if
        
        ! Copy to other climate variables
        tmp(iwst)%start_yr = pcp(iwst)%start_yr
        tmp(iwst)%end_yr = pcp(iwst)%end_yr
        tmp(iwst)%start_day = pcp(iwst)%start_day
        tmp(iwst)%end_day = pcp(iwst)%end_day
        tmp(iwst)%yrs_start = pcp(iwst)%yrs_start
        
        slr(iwst)%start_yr = pcp(iwst)%start_yr
        slr(iwst)%end_yr = pcp(iwst)%end_yr
        slr(iwst)%start_day = pcp(iwst)%start_day
        slr(iwst)%end_day = pcp(iwst)%end_day
        slr(iwst)%yrs_start = pcp(iwst)%yrs_start
        
        hmd(iwst)%start_yr = pcp(iwst)%start_yr
        hmd(iwst)%end_yr = pcp(iwst)%end_yr
        hmd(iwst)%start_day = pcp(iwst)%start_day
        hmd(iwst)%end_day = pcp(iwst)%end_day
        hmd(iwst)%yrs_start = pcp(iwst)%yrs_start
        
        wnd(iwst)%start_yr = pcp(iwst)%start_yr
        wnd(iwst)%end_yr = pcp(iwst)%end_yr
        wnd(iwst)%start_day = pcp(iwst)%start_day
        wnd(iwst)%end_day = pcp(iwst)%end_day
        wnd(iwst)%yrs_start = pcp(iwst)%yrs_start
        
        ! Allocate time series arrays
        allocate (pcp(iwst)%ts(366, pcp(iwst)%nbyr), source = 0.)
        allocate (tmp(iwst)%ts(366, tmp(iwst)%nbyr), source = 0.)   ! ts for TMAX (maximum temperature)
        allocate (tmp(iwst)%ts2(366, tmp(iwst)%nbyr), source = 0.) ! ts2 for TMIN (minimum temperature)
        allocate (slr(iwst)%ts(366, slr(iwst)%nbyr), source = 0.)
        allocate (hmd(iwst)%ts(366, hmd(iwst)%nbyr), source = 0.)
        allocate (wnd(iwst)%ts(366, wnd(iwst)%nbyr), source = 0.)
        
        ! Populate the time series data from NetCDF arrays
        ! Map NetCDF time index to (day, year) indices
        ! CRITICAL: Use Julian day indexing (1-366) like traditional method
        itime = 1
        do iyear = 1, pcp(iwst)%nbyr
            ! Handle leap years properly (using traditional method logic)
            ! Calculate actual year for this data year
            actual_year = pcp(iwst)%start_yr + iyear - 1
            days_in_year_loop = 365
            
            ! Check for leap year (same logic as time_control.f90)
            if (Mod(actual_year,4) == 0) then
              if (Mod(actual_year,100) == 0) then
                if (Mod(actual_year,400) == 0) then
                  days_in_year_loop = 366  ! Leap year
                else
                  days_in_year_loop = 365  ! Century year, not leap
                end if
              else
                days_in_year_loop = 366  ! Leap year
              end if
            else
              days_in_year_loop = 365  ! Not a leap year
            endif
            
            do iday = 1, days_in_year_loop  ! Use correct number of days for the year
                if (itime <= ntime) then
                    ! IMPORTANT: Use Julian day (iday) directly, not converted day
                    ! Traditional method: pcp(i)%ts(istep, iyrs) where istep is Julian day
                    ! NetCDF method: pcp(iwst)%ts(iday, iyear) where iday is Julian day
                    
                    ! Precipitation - use correct NetCDF indexing (lon, lat, time)
                    ! Apply station scaling factor (like traditional method)
                    if (allocated(pcp_data)) then
                        pcp(iwst)%ts(iday, iyear) = pcp_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%pcp_factor
                        if (pcp(iwst)%ts(iday, iyear) < 0.0) pcp(iwst)%ts(iday, iyear) = 0.0
                    endif
                    
                    ! Temperature - populate both ts (TMAX) and ts2 (TMIN)
                    ! Apply station scaling factors (like traditional method)
                    if (allocated(tmax_data)) then
                        tmp(iwst)%ts(iday, iyear) = tmax_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%tmax_factor
                    endif
                    
                    if (allocated(tmin_data)) then
                        tmp(iwst)%ts2(iday, iyear) = tmin_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%tmin_factor
                    endif
                    
                    ! Solar radiation
                    if (allocated(slr_data)) then
                        slr(iwst)%ts(iday, iyear) = slr_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%slr_factor
                        if (slr(iwst)%ts(iday, iyear) < 0.0) slr(iwst)%ts(iday, iyear) = 0.0
                    endif
                    
                    ! Humidity
                    if (allocated(hmd_data)) then
                        hmd(iwst)%ts(iday, iyear) = hmd_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%hmd_factor
                        if (hmd(iwst)%ts(iday, iyear) < 0.0) hmd(iwst)%ts(iday, iyear) = 0.0
                        if (hmd(iwst)%ts(iday, iyear) > 1.0) hmd(iwst)%ts(iday, iyear) = 1.0
                    endif
                    
                    ! Wind speed
                    if (allocated(wnd_data)) then
                        wnd(iwst)%ts(iday, iyear) = wnd_data(target_lon_idx, target_lat_idx, itime) * wst(iwst)%wnd_factor
                        if (wnd(iwst)%ts(iday, iyear) < 0.0) wnd(iwst)%ts(iday, iyear) = 0.0
                    endif
                    
                    itime = itime + 1
                else
                    ! No more NetCDF data available, fill with zeros or last value
                    exit
                endif
            end do
        end do
        
    end do
    
    
    ! close netcdf file
    status = nf_close(ncid)
    
    ! clean up allocated arrays
    deallocate(time_vals, lat_vals, lon_vals)
    deallocate(pcp_data, tmin_data, tmax_data, slr_data, hmd_data, wnd_data)
    
    write(*,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"
    write(9003,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"

    return

contains

    ! Helper subroutine to parse time units string
    subroutine parse_time_units(units_str, ref_yr, ref_mo, ref_dy)
        character(len=*), intent(in) :: units_str
        integer, intent(out) :: ref_yr, ref_mo, ref_dy
        
        integer :: pos1, pos2, ios
        character(len=20) :: date_part
        
        ! Initialize defaults
        ref_yr = 1970
        ref_mo = 1
        ref_dy = 1
        
        ! Find "since" keyword
        pos1 = index(units_str, "since")
        if (pos1 > 0) then
            pos1 = pos1 + 5  ! Move past "since"
            
            ! Skip whitespace
            do while (pos1 <= len(units_str) .and. units_str(pos1:pos1) == ' ')
                pos1 = pos1 + 1
            end do
            
            ! Find end of date part (before time if present)
            pos2 = index(units_str(pos1:), ' ')
            if (pos2 == 0) then
                pos2 = len(units_str) + 1
            else
                pos2 = pos1 + pos2 - 1
            endif
            
            date_part = units_str(pos1:pos2-1)
            
            ! Parse YYYY-MM-DD format
            read(date_part(1:4), *, iostat=ios) ref_yr
            if (ios == 0 .and. len_trim(date_part) >= 7) then
                read(date_part(6:7), *, iostat=ios) ref_mo
            endif
            if (ios == 0 .and. len_trim(date_part) >= 10) then
                read(date_part(9:10), *, iostat=ios) ref_dy
            endif
        endif
        
        ! Validate parsed values
        if (ref_yr < 1000 .or. ref_yr > 5000) ref_yr = 1970
        if (ref_mo < 1 .or. ref_mo > 12) ref_mo = 1
        if (ref_dy < 1 .or. ref_dy > 31) ref_dy = 1
        
    end subroutine parse_time_units
    
    ! Helper subroutine to add days to a date
    subroutine add_days_to_date(start_yr, start_mo, start_dy, days_to_add, end_yr, end_mo, end_dy)
        integer, intent(in) :: start_yr, start_mo, start_dy, days_to_add
        integer, intent(out) :: end_yr, end_mo, end_dy
        
        integer :: days_left, days_in_month
        integer :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        end_yr = start_yr
        end_mo = start_mo
        end_dy = start_dy
        days_left = days_to_add
        
        ! Add days
        do while (days_left > 0)
            ! Check for leap year and adjust February
            if (end_mo == 2) then
                if (is_leap_year(end_yr)) then
                    days_in_month = 29
                else
                    days_in_month = 28
                endif
            else
                days_in_month = month_days(end_mo)
            endif
            
            if (end_dy + days_left <= days_in_month) then
                ! Remaining days fit in current month
                end_dy = end_dy + days_left
                days_left = 0
            else
                ! Move to next month
                days_left = days_left - (days_in_month - end_dy + 1)
                end_dy = 1
                end_mo = end_mo + 1
                if (end_mo > 12) then
                    end_mo = 1
                    end_yr = end_yr + 1
                endif
            endif
        end do
        
    end subroutine add_days_to_date
    
    ! Helper function to check if year is leap year
    logical function is_leap_year(check_year)
        integer, intent(in) :: check_year
        
        is_leap_year = .false.
        if (mod(check_year, 4) == 0) then
            if (mod(check_year, 100) == 0) then
                if (mod(check_year, 400) == 0) then
                    is_leap_year = .true.
                endif
            else
                is_leap_year = .true.
            endif
        endif
        
    end function is_leap_year

end subroutine cli_ncdf_meas