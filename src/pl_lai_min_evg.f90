      function pl_lai_min_evg (j, jpl) result (lai_min)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this function returns the winter lai floor of an evergreen perennial (gro_trig = temp_gro_evg):
!!    lai_min from plants.plt, limited by the age-dependent maximum lai used in pl_leaf_gro so
!!    a young stand is never given more winter canopy than it can grow

      use plant_data_module
      use plant_module

      implicit none

      integer, intent (in) :: j         !none          |HRU number
      integer, intent (in) :: jpl       !none          |plant number in the community
      real :: lai_min                   !m**2/m**2     |winter lai floor
      integer :: idp                    !none          |plant number from plants.plt
      real :: rto_lin                   !none          |ratio of current years of growth:years to maturity
      real :: laimax                    !m**2/m**2     |age-dependent maximum lai

      idp = pcom(j)%plcur(jpl)%idplt
      rto_lin = float(max(1, pcom(j)%plcur(jpl)%curyr_mat)) / float(pldb(idp)%mat_yrs)
      laimax = pcom(j)%plcur(jpl)%lai_pot * 10. ** (alog10(rto_lin) * pldb(idp)%laixco_tree)
      laimax = min(laimax, pcom(j)%plcur(jpl)%lai_pot)
      lai_min = min(pldb(idp)%alai_min, laimax)
      lai_min = max(0., lai_min)

      return
      end function pl_lai_min_evg
