      subroutine soil_bd_update (j, ly_in, chg_typ, chg_val, src)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine changes the bulk density of one soil layer (ly_in > 0) or of every layer
!!    (ly_in = 0) of HRU j during a run and recomputes the properties derived from it
!!    (wp, up, por, ul, fc, hk, wpmm, vwt, volcr, conv_wt and the profile sums), keeping the
!!    stored layer water (st, mm above wilting point) unchanged, so soil%sw and the water balance
!!    are continuous.  awc and k are kept.  Unlike soil_awc_init, layer water is not reset to fc*ffc.
!!    if st is above the new saturation limit (ul) the excess is left in the layer, where the next
!!    day's percolation and saturation excess routines move it, and is reported in soil_bd_chg.txt

      use soil_module
      use hru_module, only : cn2
      use time_module

      implicit none

      external :: curno

      integer, intent (in) :: j                 !none          |HRU number
      integer, intent (in) :: ly_in             !none          |layer to change (0 = all layers)
      character(len=*), intent (in) :: chg_typ  !              |type of change (absval, abschg, pctchg)
      real, intent (in) :: chg_val              !              |amount of change
      character(len=*), intent (in) :: src      !              |calling operation (for the report)
      real :: chg_par                           !              |new parameter value
      integer :: ly = 0                         !none          |soil layer counter
      integer :: ly1 = 0                        !none          |first layer to change
      integer :: ly2 = 0                        !none          |last layer to change
      integer :: nly = 0                        !none          |number of soil layers
      integer, save :: iunit = 0                !none          |report file unit
      real :: sumpor = 0.                       !mm            |porosity of profile
      real :: drpor = 0.                        !              |drainable porosity
      real :: bd_old = 0.                       !Mg/m**3       |bulk density before the change
      real :: wp_old = 0.                       !mm/mm         |wilting point before the change
      real :: fc_old = 0.                       !mm            |field capacity (above wp) before the change
      real :: ul_old = 0.                       !mm            |saturation (above wp) before the change
      real :: wpmm_old = 0.                     !mm            |wilting point water before the change
      real :: excess = 0.                       !mm            |st above the new ul

      nly = soil(j)%nly
      if (ly_in > nly) return
      if (ly_in <= 0) then
        ly1 = 1
        ly2 = nly
      else
        ly1 = ly_in
        ly2 = ly_in
      end if

      if (iunit == 0) then
        open (newunit=iunit, file="soil_bd_chg.txt", recl=800)
        write (iunit,*) "soil bulk density changes during the run (st kept unchanged)"
        write (iunit,'(a6,a5,a5,a6,a9,a4,a12,15a12)') "jday", "mon", "day", "yr", "hru", "ly", "src",  &
          "bd_old", "bd_new", "wp_old", "wp_new", "wpmm_old", "wpmm_new", "fc_old", "fc_new",          &
          "ul_old", "ul_new", "st", "st_excess", "sw", "conv_wt", "por_new"
        write (iunit,'(a6,a5,a5,a6,a9,a4,a12,15a12)') "", "", "", "", "", "", "",                      &
          "Mg/m3", "Mg/m3", "mm/mm", "mm/mm", "mm", "mm", "mm", "mm",                                 &
          "mm", "mm", "mm", "mm", "mm", "kg/ha/ppm", "mm/mm"
      end if

      do ly = ly1, ly2
        bd_old = soil(j)%phys(ly)%bd
        wp_old = soil(j)%phys(ly)%wp
        fc_old = soil(j)%phys(ly)%fc
        ul_old = soil(j)%phys(ly)%ul
        wpmm_old = soil(j)%phys(ly)%wpmm

        soil(j)%phys(ly)%bd = chg_par (soil(j)%phys(ly)%bd, chg_typ, chg_val, 0.9, 2.5)

        !! water contents at -1.5 MPa and -0.033 MPa and porosity (as soil_awc_init)
        if (soil(j)%phys(ly)%awc <= 1.e-6) soil(j)%phys(ly)%awc = .005
        if (soil(j)%phys(ly)%awc >= .8) soil(j)%phys(ly)%awc = .8
        soil(j)%phys(ly)%wp = 0.4 * soil(j)%phys(ly)%clay * soil(j)%phys(ly)%bd / 100.
        if (soil(j)%phys(ly)%wp <= 0.) soil(j)%phys(ly)%wp = .005
        soil(j)%phys(ly)%up = soil(j)%phys(ly)%wp + soil(j)%phys(ly)%awc
        soil(j)%phys(ly)%por = 1. - soil(j)%phys(ly)%bd / 2.65
        if (soil(j)%phys(ly)%up >= soil(j)%phys(ly)%por) then
          soil(j)%phys(ly)%up = soil(j)%phys(ly)%por - .05
          soil(j)%phys(ly)%wp = soil(j)%phys(ly)%up - soil(j)%phys(ly)%awc
          if (soil(j)%phys(ly)%wp <= 0.) then
            soil(j)%phys(ly)%up = soil(j)%phys(ly)%por * .75
            soil(j)%phys(ly)%wp = soil(j)%phys(ly)%por * .25
          end if
        end if
        drpor = soil(j)%phys(ly)%por - soil(j)%phys(ly)%up
        soil(j)%ly(ly)%vwt = (437.13 * drpor**2) - (95.08 * drpor) + 8.257

        !! layer capacities; st is not changed
        soil(j)%phys(ly)%ul = (soil(j)%phys(ly)%por - soil(j)%phys(ly)%wp) * soil(j)%phys(ly)%thick
        soil(j)%phys(ly)%fc = soil(j)%phys(ly)%thick * (soil(j)%phys(ly)%up - soil(j)%phys(ly)%wp)
        soil(j)%phys(ly)%hk = (soil(j)%phys(ly)%ul - soil(j)%phys(ly)%fc) / soil(j)%phys(ly)%k
        if (soil(j)%phys(ly)%hk < 1.) soil(j)%phys(ly)%hk = 1.
        soil(j)%phys(ly)%wpmm = soil(j)%phys(ly)%wp * soil(j)%phys(ly)%thick
        if (soil(j)%phys(ly)%fc > 0.) then
          soil(j)%ly(ly)%volcr = soil(j)%phys(ly)%crdep * (soil(j)%phys(ly)%fc - soil(j)%phys(ly)%st) / &
                                                                                   soil(j)%phys(ly)%fc
        end if
        soil(j)%phys(ly)%conv_wt = soil(j)%phys(ly)%bd * soil(j)%phys(ly)%thick / 100.

        excess = max (0., soil(j)%phys(ly)%st - soil(j)%phys(ly)%ul)
        write (iunit,'(i6,i5,i5,i6,i9,i4,a12,15f12.4)') time%day, time%mo, time%day_mo, time%yrc, j, ly,    &
          trim(src), bd_old, soil(j)%phys(ly)%bd, wp_old, soil(j)%phys(ly)%wp, wpmm_old,                    &
          soil(j)%phys(ly)%wpmm, fc_old, soil(j)%phys(ly)%fc, ul_old, soil(j)%phys(ly)%ul,                  &
          soil(j)%phys(ly)%st, excess, soil(j)%sw, soil(j)%phys(ly)%conv_wt, soil(j)%phys(ly)%por
      end do

      !! profile sums (sw is unchanged because st is unchanged)
      sumpor = 0.
      soil(j)%sumfc = 0.
      soil(j)%sumul = 0.
      soil(j)%sumwp = 0.
      do ly = 1, nly
        sumpor = sumpor + soil(j)%phys(ly)%por * soil(j)%phys(ly)%thick
        soil(j)%sumfc = soil(j)%sumfc + soil(j)%phys(ly)%fc
        soil(j)%sumul = soil(j)%sumul + soil(j)%phys(ly)%ul
        soil(j)%sumwp = soil(j)%sumwp + soil(j)%phys(ly)%wpmm
      end do
      soil(j)%avpor = sumpor / soil(j)%phys(nly)%d
      soil(j)%avbd = 2.65 * (1. - soil(j)%avpor)

      !! curve number coefficients depend on sumul and sumfc
      call curno (cn2(j), j)

      return
      end subroutine soil_bd_update
