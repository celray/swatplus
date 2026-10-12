#!/usr/bin/env python3
"""Check the NetCDF reader's array population with gfortran or ifx.

Run: python3 test/test_netcdf_missing.py [--compiler gfortran] [--source PATH]
The production routine is compiled unchanged with small climate arrays.
"""
import argparse
from pathlib import Path
import subprocess
import tempfile

PREFIX = """
module climate_check
use iso_c_binding, only: c_float
implicit none
type series
 real :: ts(366,1)=0., ts2(366,1)=0.
 integer :: nbyr=1, start_yr=2001
end type
type station
 real :: pcp_factor=1.,tmax_factor=1.,tmin_factor=1.,slr_factor=1.,hmd_factor=1.,wnd_factor=1.
end type
type(series) :: pcp(1),tmp(1),slr(1),hmd(1),wnd(1)
type(station) :: wst(1)
real(c_float),allocatable :: pcp_data(:,:,:),tmax_data(:,:,:),tmin_data(:,:,:)
real(c_float),allocatable :: slr_data(:,:,:),hmd_data(:,:,:),wnd_data(:,:,:)
integer :: itime,itime_offset=0,iyear,actual_year,days_in_year_loop,iday
contains
logical function is_leap_year(y)
integer,intent(in)::y
is_leap_year=mod(y,4)==0.and.(mod(y,100)/=0.or.mod(y,400)==0)
end function
"""
SUFFIX = """
end module
program check
use climate_check
implicit none
real :: missing(6)=[-97.,-98.,-99.,-100.,-999.,-9999.]
real :: factors(6)=[1.,0.001,0.,-1.,2.,-99.]
real :: valid(8)=[-96.,-20.,-0.5,-0.,0.,0.5,1.,23.4567]
real :: actual(6),expected(6),value,factor
integer :: i,j,k,checks
allocate(pcp_data(1,1,1),tmax_data(1,1,1),tmin_data(1,1,1))
allocate(slr_data(1,1,1),hmd_data(1,1,1),wnd_data(1,1,1))
checks=0
do i=1,size(missing)
 do j=1,size(factors)
  value=missing(i); factor=factors(j)
  call populate(value,factor,actual)
  if(any(actual/=value)) then
   print *, 'Missing value changed:',value,'factor:',factor,'stored:',actual
   stop 1
  endif
  if(any(actual > -97.)) stop 2
  checks=checks+6
 enddo
enddo
! Compare valid values bit for bit with the previous scaling and bounds.
do i=1,size(valid)
 do j=1,size(factors)
  value=valid(i); factor=factors(j)
  expected=value*factor
  do k=1,6
   if(k==2.or.k==3) cycle
   if(expected(k)<0.) expected(k)=0.
  enddo
  if(expected(5)>1.) expected(5)=1.
  call populate(value,factor,actual)
  if(any(transfer(actual,[0,0,0,0,0,0])/=transfer(expected,[0,0,0,0,0,0]))) then
   print *, 'Valid value changed:',value,'factor:',factor
   stop 3
  endif
  checks=checks+6
 enddo
enddo
! One missing variable must leave the other daily variables unchanged.
pcp_data=1.;tmax_data=-999.;tmin_data=2.;slr_data=3.;hmd_data=0.5;wnd_data=4.
wst(1)=station()
call populate_timeseries_data(1,1,1,1)
call stored(actual)
expected=[1.,-999.,2.,3.,0.5,4.]
if(any(actual/=expected)) stop 4
checks=checks+6
print *, checks, 'array-value checks passed'
contains
subroutine stored(values)
 real,intent(out)::values(6)
 values=[pcp(1)%ts(1,1),tmp(1)%ts(1,1),tmp(1)%ts2(1,1), &
         slr(1)%ts(1,1),hmd(1)%ts(1,1),wnd(1)%ts(1,1)]
end subroutine
subroutine populate(value,factor,values)
 real,intent(in)::value,factor
 real,intent(out)::values(6)
 pcp_data=value;tmax_data=value;tmin_data=value;slr_data=value;hmd_data=value;wnd_data=value
 wst(1)=station(factor,factor,factor,factor,factor,factor)
 call populate_timeseries_data(1,1,1,1)
 call stored(values)
end subroutine
end program
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] / 'src/cli_ncdf_meas.f90')
    args = parser.parse_args()
    source = args.source.read_text()
    start = source.index('    subroutine populate_timeseries_data(')
    end_marker = '    end subroutine populate_timeseries_data'
    end = source.index(end_marker, start) + len(end_marker)
    flags = ['-O2', '-ffree-line-length-none', '-fcheck=all']
    if Path(args.compiler).name == 'ifx':
        flags = ['-O2', '-free', '-check', 'all']
    with tempfile.TemporaryDirectory(prefix='netcdf-missing-') as folder:
        root = Path(folder)
        (root / 'check.f90').write_text(PREFIX + source[start:end] + SUFFIX)
        subprocess.run([args.compiler, *flags, 'check.f90', '-o', 'check'], cwd=root, check=True)
        subprocess.run([str(root / 'check')], cwd=root, check=True)


if __name__ == '__main__':
    main()
