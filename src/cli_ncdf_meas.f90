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
    integer :: days_since_1970, year, month, day, days_in_months(12)
    
    ! Loop counters and diagnostics
    integer :: istat, itime, iyear, iday, i, iwst
    integer :: start_year, start_day, current_year, days_in_year
    integer :: actual_year  ! For leap year calculation
    integer :: days_in_year_loop ! For proper leap year counting
    integer :: ndims, dimids(NF_MAX_VAR_DIMS), natts, xtype
    character(len=NF_MAX_NAME) :: var_name, dim_name
    integer :: dim_len
    
    version_string = nf_inq_libvers()
    write (*,*) "reading data using netcdf version", trim(version_string)
    write (9003,*) "reading data using netcdf version", trim(version_string)

    ! Get NetCDF filename - use test data file
    ncdf_file = "save.nc4"
    
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
    
    
    ! Find closest grid point
    min_dist = 999999.0
    target_lat_idx = 1
    target_lon_idx = 1
    
    do ilat = 1, nlat
        do ilon = 1, nlon
            ! Simple Euclidean distance (for small areas this is adequate)
            dist = sqrt((lat_vals(ilat) - target_lat)**2 + (lon_vals(ilon) - target_lon)**2)
            if (dist < min_dist) then
                min_dist = dist
                target_lat_idx = ilat
                target_lon_idx = ilon
            endif
        end do
    end do
    
100 continue
    
    ! Read time data to get the date for first time step
    status = nf_inq_varid(ncid, "time", time_varid)
    if (status == NF_NOERR) then
        status = nf_get_var_real(ncid, time_varid, time_vals)
        if (status == NF_NOERR) then
            if (ntime >= 1) then
                ! Time is in "days since 1970-01-01 00:00:00.0"
                ! Convert to readable date
                days_since_1970 = int(time_vals(1))
                
                ! Simple date calculation (approximate)
                year = 1970 + days_since_1970 / 365
                ! For simplicity, just show the raw days value and approximate year
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
        if (Mod(actual_year, 4) == 0) then
            pcp(iwst)%end_day = 366  ! Last day of leap year
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
            ! Handle leap years properly (like traditional method)
            ! Calculate actual year for this data year
            actual_year = pcp(iwst)%start_yr + iyear - 1
            days_in_year_loop = 365
            
            ! Check for leap year (same logic as traditional method)
            if (Mod(actual_year, 4) == 0) then
                days_in_year_loop = 366
            else
                days_in_year_loop = 365
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
    
    write(*,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"
    write(9003,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"

    
    ! close netcdf file
    status = nf_close(ncid)
    
    ! clean up allocated arrays
    deallocate(time_vals, lat_vals, lon_vals)
    deallocate(pcp_data, tmin_data, tmax_data, slr_data, hmd_data, wnd_data)

    write (*,*) "netcdf weather data reading completed successfully!"
    write (9003,*) "netcdf weather data reading completed successfully!"

    return

end subroutine cli_ncdf_meas