      subroutine soil_bd_update (j, ly_str, chg_typ, chg_val, src, op_name)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine changes the bulk density of one soil layer or of every layer
!!    ("all", "null" or "0") of HRU j during a run and recomputes the properties derived from it
!!    (wp, up, por, ul, fc, hk, wpmm, vwt, volcr, conv_wt and the profile sums), keeping the
!!    stored layer water (st, mm above wilting point) unchanged, so soil%sw and the water balance
!!    are continuous.  awc and k are kept.  Unlike soil_awc_init, layer water is not reset to fc*ffc.
!!    if st is above the new saturation limit (ul) the excess is left in the layer, where the
!!    percolation and saturation excess routines move it (the same day when the change is made by a
!!    decision table or schedule), and is reported in soil_bd_chg.txt
!!    water below the wilting point is not a stored state: keeping st changes it by wpmm_old - wpmm_new,
!!    reported as wpmm_chg in soil_bd_chg.txt
!!    bd is limited to 0.9 - 2.0 Mg/m3: 2.0 is the clamp in soil_phys_init and keeps the soil
!!    temperature water factor in stmp_solt (0.356 - 0.144 * avbd) positive; 0.9 is the cal_parms.cal minimum
!!    the run stops on a change type other than absval, abschg or pctchg, or on a layer that is not
!!    "all", "null", "0" or a layer number from 1 to the number of layers

      use soil_module
      use hru_module, only : cn2
      use time_module

      implicit none

      external :: curno

      integer, intent (in) :: j                 !none          |HRU number
      character(len=*), intent (in) :: ly_str   !              |layer to change ("all", "null" or "0" = all layers)
      character(len=*), intent (in) :: chg_typ  !              |type of change (absval, abschg, pctchg)
      real, intent (in) :: chg_val              !              |amount of change
      character(len=*), intent (in) :: src      !              |calling operation (bd_update or bdup)
      character(len=*), intent (in) :: op_name  !              |name of the action or operation (for messages)
      integer :: ly_in = 0                      !none          |layer to change (0 = all layers)
      integer :: ios = 0                        !none          |read status of the layer
      character(len=16) :: src_file = ""        !              |input file of the calling operation
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
      if (src == "bdup") then
        src_file = "management.sch"
      else
        src_file = "lum.dtl"
      end if

      !! change type
      if (chg_typ /= "absval" .and. chg_typ /= "abschg" .and. chg_typ /= "pctchg") then
        write (*,'(5a,i0,3a)') "Error: ", trim(src), " ", trim(op_name), " for HRU ", j, " in ", trim(src_file),  &
          " has change type '" // trim(chg_typ) // "'; use absval, abschg or pctchg"
        print *
        error stop
      end if

      !! layer
      if (ly_str == "all" .or. ly_str == "null" .or. ly_str == "0") then
        ly_in = 0
      else
        read (ly_str,*,iostat=ios) ly_in
        if (ios /= 0 .or. ly_in < 1 .or. ly_in > nly) then
          write (*,'(5a,i0,3a,i0,a)') "Error: ", trim(src), " ", trim(op_name), " for HRU ", j, " in ", trim(src_file), &
            " has layer '" // trim(ly_str) // "'; use all, null, 0 or a layer from 1 to ", nly, " for this soil"
          print *
          error stop
        end if
      end if

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
        write (iunit,'(a6,a5,a5,a6,a9,a4,a12,16a12)') "jday", "mon", "day", "yr", "hru", "ly", "src",  &
          "bd_old", "bd_new", "wp_old", "wp_new", "wpmm_old", "wpmm_new", "fc_old", "fc_new",          &
          "ul_old", "ul_new", "st", "st_excess", "sw", "conv_wt", "por_new", "wpmm_chg"
        write (iunit,'(a6,a5,a5,a6,a9,a4,a12,16a12)') "", "", "", "", "", "", "",                      &
          "Mg/m3", "Mg/m3", "mm/mm", "mm/mm", "mm", "mm", "mm", "mm",                                 &
          "mm", "mm", "mm", "mm", "mm", "kg/ha/ppm", "mm/mm", "mm"
      end if

      do ly = ly1, ly2
        bd_old = soil(j)%phys(ly)%bd
        wp_old = soil(j)%phys(ly)%wp
        fc_old = soil(j)%phys(ly)%fc
        ul_old = soil(j)%phys(ly)%ul
        wpmm_old = soil(j)%phys(ly)%wpmm

        soil(j)%phys(ly)%bd = chg_par (soil(j)%phys(ly)%bd, chg_typ, chg_val, 0.9, 2.0)

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
        write (iunit,'(i6,i5,i5,i6,i9,i4,a12,16f12.4)') time%day, time%mo, time%day_mo, time%yrc, j, ly,    &
          trim(src), bd_old, soil(j)%phys(ly)%bd, wp_old, soil(j)%phys(ly)%wp, wpmm_old,                    &
          soil(j)%phys(ly)%wpmm, fc_old, soil(j)%phys(ly)%fc, ul_old, soil(j)%phys(ly)%ul,                  &
          soil(j)%phys(ly)%st, excess, soil(j)%sw, soil(j)%phys(ly)%conv_wt, soil(j)%phys(ly)%por,          &
          wpmm_old - soil(j)%phys(ly)%wpmm
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
