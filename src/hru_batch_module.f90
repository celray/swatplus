      module hru_batch_module

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    runs groups (batches) of hru commands in parallel with OpenMP and gives
!!    the same results as the serial command loop
!!
!!    a batch is a run of consecutive hru commands in which no hru receives
!!    flow from another hru of the same run. the hrus of a batch are given
!!    slots 1, 2, ... in command order. code that touches state another hru
!!    of the batch can change (water sources, files, basin totals, decision
!!    tables shared between hrus) calls ordered_gate first. the gate waits
!!    until every earlier slot has finished, so that code sees the same state
!!    as in the serial loop.
!!
!!    without OpenMP the batches run one hru after the other and the gate
!!    returns at once.

      use hydrograph_module, only : hyd_output
      use gwflow_defer_module, only : gwdef_init, gwdef_apply
      use conditional_module, only : decision_table
      use mgt_operations_module, only : management_ops
      use organic_mineral_mass_module, only : organic_mass

      implicit none

      integer :: nbatch = 0                                 !none       |number of batches
      integer :: batch_min = 2                              !none       |fewest hrus in a batch
      integer, dimension (:), allocatable :: batch_beg      !none       |first position of batch in batch_ob
      integer, dimension (:), allocatable :: batch_end      !none       |last position of batch in batch_ob
      integer, dimension (:), allocatable :: batch_ob       !none       |object numbers of all batches in command order
      integer, dimension (:), allocatable :: ob_batch       !none       |batch that starts at this object (0 = none)
      logical :: batch_on = .false.                         !none       |true if batches are used

      !! ordered gate
      integer :: islot = 0                                  !none       |slot of the hru this thread is running
      integer :: slots_done = 0                             !none       |every slot up to this one has finished
      logical, dimension (:), allocatable :: slot_done      !none       |slot has finished
      type (hyd_output), dimension (:), allocatable :: slot_ht1   !   |deposition hydrograph (ht1) of each slot
      !$omp threadprivate (islot)

      !! decision tables. a table that can touch state of other hrus (dtbl_coupled) runs
      !! on the shared table after the gate; the others run on a copy owned by the thread
      logical, dimension (:), allocatable :: lum_coupled    !none       |dtbl_lum table is coupled
      logical, dimension (:), allocatable :: res_coupled    !none       |dtbl_res table is coupled
      type (decision_table), dimension (:), allocatable, target :: lum_copy   !  |thread copy of dtbl_lum
      type (decision_table), dimension (:), allocatable, target :: res_copy   !  |thread copy of dtbl_res
      !$omp threadprivate (lum_copy, res_copy)
      logical :: graze_gate = .false.                       !none       |more than one grazing operation

      !! milestones: points in hru_control after which a slot no longer changes one of
      !! the leftover values. a slot that needs the value of an earlier slot waits with
      !! ordered_wait until every earlier slot has passed the milestone (or finished),
      !! which is sooner than ordered_gate
      integer, parameter :: nmile = 2                       !none       |milestones
      integer, parameter :: mile_enratio = 1                !none       |enratio is final (after pest_enrsb)
      integer, parameter :: mile_plgrow = 2                 !none       |delg_p is final (after pl_grow)
      integer, dimension (nmile), parameter :: mile_real = (/ 6, 8 /)   !   |leftover real of each milestone (leftover_real)
      logical, dimension (:,:), allocatable :: slot_mile    !none       |slot passed the milestone
      integer, dimension (nmile) :: miles_done = 0          !none       |every slot up to this one passed it
      integer, dimension (nmile) :: mile_seen = 0           !none       |slots searched for the last setter
      integer, dimension (:,:), allocatable :: mile_setter  !none       |last slot up to this slot that set the value

      !! statistics (SWATPLUS_BATCH_STATS=1, printed by batch_stats_print)
      logical :: stats_on = .false.                         !none       |collect statistics
      integer :: nsite = 0                                  !none       |gate call sites seen
      character(len=16), dimension (32) :: site_name = ""   !           |gate call site
      integer (kind=8), dimension (32) :: site_calls = 0    !           |gate calls
      integer (kind=8), dimension (32) :: site_waits = 0    !           |gate calls that waited
      integer (kind=8), dimension (32) :: site_ticks = 0    !           |clock ticks spent in the gate
      integer (kind=8) :: run_ticks = 0                     !           |clock ticks in the parallel part of batches
      integer (kind=8) :: post_ticks = 0                    !           |clock ticks in the ordered part of batches

      !! leftover values: values an hru can read before it writes them. the serial loop
      !! then uses what the closest earlier hru (or call) left, so each slot starts with
      !! the marker unset, publishes what it wrote in slot_finish, and ordered_gate fills
      !! in the rest from the earlier slots (slot 0 is the value before the batch).
      !! the lo_ variables are locals of the routines named below that keep their value
      !! between calls in the serial code; the routines use them under the local name
      real, parameter :: unset = -huge(1.)                  !none       |marker for not written (real)
      integer, parameter :: iunset = -huge(1)               !none       |marker for not written (integer)
      real :: lo_delg_p = 0.                                !           |pl_nut_demand delg_p
      real :: lo_c_microb_sed = 0.                          !           |nut_orgnc2 c_microb_sed
      real :: lo_gsi_adj = 0.                               !           |et_pot gsi_adj
      real :: lo_yield = 0.                                 !           |actions yield
      integer :: lo_j = 0                                   !           |actions j
      integer :: lo_irrop = 0                               !           |actions irrop
      integer :: lo_idtill = 0                              !           |actions idtill
      integer :: lo_istr = 0                                !           |actions istr
      integer :: lo_istr1 = 0                               !           |actions istr1
      !$omp threadprivate (lo_delg_p, lo_c_microb_sed, lo_gsi_adj, lo_yield, lo_j, lo_irrop, &
      !$omp& lo_idtill, lo_istr, lo_istr1)
      integer, parameter :: nreal = 11                      !none       |leftover reals (leftover_get)
      integer, parameter :: nint = 5                        !none       |leftover integers (leftover_get)
      integer, parameter :: narr = 3                        !none       |leftover plant arrays (un2, up2, uapd)
      type (management_ops), dimension (:), allocatable :: slot_mgt
      type (organic_mass), dimension (:), allocatable :: slot_mass_up
      real, dimension (:,:), allocatable :: slot_real       !           |leftover reals (value, slot)
      integer, dimension (:,:), allocatable :: slot_int     !           |leftover integers (value, slot)
      real, dimension (:,:,:), allocatable :: slot_arr      !           |leftover arrays (array, plant, slot)
      logical, dimension (:), allocatable :: slot_mgt_set
      logical, dimension (:), allocatable :: slot_mass_up_set
      logical, dimension (:,:), allocatable :: slot_real_set
      logical, dimension (:,:), allocatable :: slot_int_set
      logical, dimension (:,:,:), allocatable :: slot_arr_set

      !! basin sums the hrus add to, kept per hru and added in command order by hru_sums
      real, dimension (:), allocatable :: wet_dep_d         !tons       |wetland deposition of the day
      logical, dimension (:), allocatable :: wet_dep_on     !none       |wet_dep_d was set today

      contains

!!    build the batches from the command sequence
      subroutine batch_init

      use hydrograph_module, only : ob, sp_ob, sp_ob1
      use maximum_data_module, only : db_mx
      use conditional_module, only : dtbl_lum, dtbl_res
      use hru_module, only : epmax, htfac, par, translt, uapd, un2, uno3d, up2, ubnrunoff, ubntss, wnan

      integer :: icmd                 !none       |command (object) number
      integer :: iob                  !none       |incoming object
      integer :: in                   !none       |counter
      integer :: npos                 !none       |positions used in batch_ob
      integer :: nrun                 !none       |hrus in the current run
      integer :: run_first            !none       |position of the first hru of the current run
      integer :: mx                   !none       |most slots in a batch
      integer :: npl                  !none       |size of the plant scratch arrays
      integer :: nts                  !none       |size of the time step scratch arrays
      integer :: nwn                  !none       |size of wnan
      integer, dimension (:), allocatable :: in_run   !none |run number each object belongs to (0 = none)
      character(len=16) :: env        !           |value of SWATPLUS_BATCH
      integer :: env_len              !           |length of env
      integer :: env_stat             !           |status of env

      allocate (wet_dep_d(sp_ob%hru), source = 0.)
      allocate (wet_dep_on(sp_ob%hru), source = .false.)
      call gwdef_init

      env = ""
      call get_environment_variable ("SWATPLUS_BATCH_STATS", env, env_len, env_stat)
      stats_on = env == "1" .or. env == "on"

      batch_on = .false.
      nbatch = 0
      env = ""
      call get_environment_variable ("SWATPLUS_BATCH", env, env_len, env_stat)
      if (env == "off" .or. env == "0") return
      if (sp_ob%hru < batch_min) return

      allocate (ob_batch(sp_ob%objs), source = 0)
      allocate (in_run(sp_ob%objs), source = 0)
      allocate (batch_ob(sp_ob%hru), source = 0)
      allocate (batch_beg(sp_ob%hru), source = 0)
      allocate (batch_end(sp_ob%hru), source = 0)

      npos = 0
      nrun = 0
      run_first = 1
      icmd = sp_ob1%objs
      do while (icmd /= 0)
        if (ob(icmd)%typ == "hru") then
          !! end the run if this hru receives from an hru of the run
          do in = 1, ob(icmd)%rcv_tot
            iob = ob(icmd)%obj_in(in)
            if (in_run(iob) == nbatch + 1 .and. nrun > 0) then
              call batch_close (run_first, npos, nrun)
              exit
            end if
          end do
          if (nrun == 0) run_first = npos + 1
          npos = npos + 1
          nrun = nrun + 1
          batch_ob(npos) = icmd
          in_run(icmd) = nbatch + 1
        else
          if (nrun > 0) call batch_close (run_first, npos, nrun)
        end if
        icmd = ob(icmd)%cmd_next
      end do
      if (nrun > 0) call batch_close (run_first, npos, nrun)

      mx = 0
      do in = 1, nbatch
        mx = max (mx, batch_end(in) - batch_beg(in) + 1)
      end do
      allocate (slot_done(mx), source = .false.)
      allocate (slot_mile(nmile,mx), source = .false.)
      allocate (mile_setter(nmile,0:mx), source = 0)
      allocate (slot_ht1(mx))
      allocate (slot_mgt(0:mx))
      allocate (slot_mass_up(0:mx))
      allocate (slot_real(nreal,0:mx), source = 0.)
      allocate (slot_int(nint,0:mx), source = 0)
      allocate (slot_arr(narr,size(uapd),0:mx), source = 0.)
      allocate (slot_mgt_set(0:mx), source = .false.)
      allocate (slot_mass_up_set(0:mx), source = .false.)
      allocate (slot_real_set(nreal,0:mx), source = .false.)
      allocate (slot_int_set(nint,0:mx), source = .false.)
      allocate (slot_arr_set(narr,size(uapd),0:mx), source = .false.)

      allocate (lum_coupled(db_mx%dtbl_lum), source = .true.)
      do in = 1, db_mx%dtbl_lum
        lum_coupled(in) = dtbl_coupled (dtbl_lum(in))
      end do
      allocate (res_coupled(db_mx%dtbl_res), source = .true.)
      do in = 1, db_mx%dtbl_res
        res_coupled(in) = dtbl_coupled (dtbl_res(in))
      end do
      graze_gate = db_mx%grazeop_db > 1

      !! every thread allocates its copy of the scratch arrays (with ifx 2026.0, copyin
      !! of an allocatable array leaves the threads on the main thread's data)
      npl = size (uapd)
      nts = size (ubnrunoff)
      nwn = size (wnan)
      !$omp parallel
      if (.not. allocated (uapd)) then
        allocate (epmax(npl), source = 0.)
        allocate (htfac(npl), source = 0.)
        allocate (par(npl), source = 0.)
        allocate (translt(npl), source = 0.)
        allocate (uapd(npl), source = 0.)
        allocate (un2(npl), source = 0.)
        allocate (uno3d(npl), source = 0.)
        allocate (up2(npl), source = 0.)
        allocate (ubnrunoff(nts), source = 0.)
        allocate (ubntss(nts), source = 0.)
        allocate (wnan(nwn), source = 0.)
      end if
      !$omp end parallel
      batch_on = nbatch > 0

      deallocate (in_run)
      return
      end subroutine batch_init

!!    close the current run; keep it as a batch if it is long enough
      subroutine batch_close (run_first, npos, nrun)

      integer, intent (in) :: run_first       !none   |position of the first hru of the run
      integer, intent (inout) :: npos         !none   |positions used in batch_ob
      integer, intent (inout) :: nrun         !none   |hrus in the run

      if (nrun >= batch_min) then
        nbatch = nbatch + 1
        batch_beg(nbatch) = run_first
        batch_end(nbatch) = npos
        ob_batch(batch_ob(run_first)) = nbatch
      else
        !! too short - these hrus stay in the serial loop
        npos = run_first - 1
      end if
      nrun = 0
      return
      end subroutine batch_close

!!    wait until every slot before this thread's slot has finished
      subroutine ordered_gate (site)

      character(len=*), intent (in), optional :: site     !   |name of the call site (statistics)
      integer :: done                 !none       |slots finished in order
      integer :: nwait                !none       |times the wait loop ran
      integer (kind=8) :: t0          !           |clock at the start of the wait
      integer (kind=8) :: t1          !           |clock at the end of the wait

      if (islot == 0) return
      nwait = 0
      if (stats_on) call system_clock (t0)
      do
        !$omp atomic read
        done = slots_done
        if (done >= islot - 1) exit
        nwait = nwait + 1
        !$omp taskyield
      end do
      !$omp flush
      if (stats_on) then
        call system_clock (t1)
        if (present (site)) then
          call stats_gate (site, nwait > 0, t1 - t0)
        else
          call stats_gate ("unnamed", nwait > 0, t1 - t0)
        end if
      end if
      call leftover_load (islot)
      return
      end subroutine ordered_gate

!!    add a gate call to the statistics
      subroutine stats_gate (site, waited, ticks)

      character(len=*), intent (in) :: site       !   |name of the call site
      logical, intent (in) :: waited              !   |the slot had to wait
      integer (kind=8), intent (in) :: ticks      !   |clock ticks spent in the gate
      integer :: i                                !   |site

      !$omp critical (batch_stats)
      do i = 1, nsite
        if (site_name(i) == site) exit
      end do
      if (i > nsite .and. nsite < size(site_name)) then
        nsite = nsite + 1
        site_name(nsite) = site
        i = nsite
      end if
      if (i <= nsite) then
        site_calls(i) = site_calls(i) + 1
        if (waited) site_waits(i) = site_waits(i) + 1
        site_ticks(i) = site_ticks(i) + ticks
      end if
      !$omp end critical (batch_stats)
      return
      end subroutine stats_gate

!!    print the batch statistics to the console (SWATPLUS_BATCH_STATS=1)
      subroutine batch_stats_print

      use iso_fortran_env, only : error_unit

      integer (kind=8) :: rate        !           |clock ticks per second
      integer :: i                    !none       |site

      if (.not. stats_on) return
      call system_clock (count_rate = rate)
      write (error_unit,'(a,i0,a,i0)') " hru batches: ", nbatch, ", largest: ", size(slot_done)
      write (error_unit,'(a,f10.2,a)') " time in batches (parallel part): ", real(run_ticks) / rate, " s"
      write (error_unit,'(a,f10.2,a)') " time in batches (ordered part):  ", real(post_ticks) / rate, " s"
      write (error_unit,'(a)') " gate site            calls      waited   seconds in gate (all threads)"
      do i = 1, nsite
        write (error_unit,'(1x,a16,2i12,f12.2)') site_name(i), site_calls(i), site_waits(i), real(site_ticks(i)) / rate
      end do
      return
      end subroutine batch_stats_print

!!    mark a slot finished and move the count of slots finished in order
      subroutine slot_finish (k)

      integer, intent (in) :: k       !none       |slot
      integer :: done                 !none       |slots finished in order

      call leftover_save (k)
      !$omp flush
      !$omp critical (batch_slots)
      slot_done(k) = .true.
      done = slots_done
      do while (done < size(slot_done))
        if (.not. slot_done(done + 1)) exit
        done = done + 1
      end do
      !$omp atomic write
      slots_done = done
      call miles_move
      !$omp end critical (batch_slots)
      return
      end subroutine slot_finish

!!    the slot running on this thread has passed milestone m: publish the leftover
!!    value of the milestone and move the count of slots past it
      subroutine milestone (m)

      integer, intent (in) :: m       !none       |milestone
      integer :: i                    !none       |leftover real of the milestone

      if (islot == 0) return
      i = mile_real(m)
      slot_real(i,islot) = leftover_real (i)
      slot_real_set(i,islot) = slot_real(i,islot) /= unset
      !$omp flush
      !$omp critical (batch_slots)
      slot_mile(m,islot) = .true.
      call miles_move
      !$omp end critical (batch_slots)
      return
      end subroutine milestone

!!    move the counts of slots past each milestone (called inside critical batch_slots)
      subroutine miles_move

      integer :: m                    !none       |milestone
      integer :: done                 !none       |slots past the milestone in order

      do m = 1, nmile
        done = miles_done(m)
        do while (done < size(slot_done))
          if (.not. (slot_mile(m,done+1) .or. slot_done(done+1))) exit
          done = done + 1
        end do
        !$omp atomic write
        miles_done(m) = done
      end do
      return
      end subroutine miles_move

!!    wait until every earlier slot has passed milestone m, then give the leftover
!!    value of the milestone the value of the closest earlier slot that set it
      subroutine ordered_wait (m, site)

      integer, intent (in) :: m                           !none   |milestone
      character(len=*), intent (in) :: site               !   |name of the call site (statistics)
      integer :: done                 !none       |slots past the milestone
      integer :: nwait                !none       |times the wait loop ran
      integer :: i                    !none       |leftover real of the milestone
      integer :: kk                   !none       |earlier slot
      integer (kind=8) :: t0          !           |clock at the start of the wait
      integer (kind=8) :: t1          !           |clock at the end of the wait

      if (islot == 0) return
      nwait = 0
      if (stats_on) call system_clock (t0)
      do
        !$omp atomic read
        done = miles_done(m)
        if (done >= islot - 1) exit
        nwait = nwait + 1
        !$omp taskyield
      end do
      !$omp flush
      if (stats_on) then
        call system_clock (t1)
        call stats_gate (site, nwait > 0, t1 - t0)
      end if
      i = mile_real(m)
      if (leftover_real (i) == unset) then
        !! the last earlier slot that set the value; the search goes on from where
        !! the previous one stopped, since all slots before islot have published
        !$omp critical (batch_slots)
        do kk = mile_seen(m) + 1, islot - 1
          mile_setter(m,kk) = mile_setter(m,kk-1)
          if (slot_real_set(i,kk)) mile_setter(m,kk) = kk
        end do
        mile_seen(m) = max (mile_seen(m), islot - 1)
        kk = mile_setter(m,islot-1)
        !$omp end critical (batch_slots)
        call leftover_real_set (i, slot_real(i,kk))
      end if
      return
      end subroutine ordered_wait

!!    one leftover real of this thread
      real function leftover_real (i)

      use hru_module, only : enratio

      integer, intent (in) :: i                           !none   |leftover real

      select case (i)
      case (6)
        leftover_real = enratio
      case (8)
        leftover_real = lo_delg_p
      case default
        leftover_real = unset
      end select
      return
      end function leftover_real

!!    set one leftover real of this thread
      subroutine leftover_real_set (i, v)

      use hru_module, only : enratio

      integer, intent (in) :: i                           !none   |leftover real
      real, intent (in) :: v                              !       |value

      select case (i)
      case (6)
        enratio = v
      case (8)
        lo_delg_p = v
      end select
      return
      end subroutine leftover_real_set

!!    reset the gate before a batch runs
      subroutine slots_reset (nslot)

      integer, intent (in) :: nslot   !none       |slots in the batch

      slot_done(1:nslot) = .false.
      slot_done(nslot+1:) = .true.
      slots_done = 0
      slot_mile(:,1:nslot) = .false.
      slot_mile(:,nslot+1:) = .true.
      miles_done = 0
      mile_seen = 0
      mile_setter(:,0) = 0
      islot = 0
      return
      end subroutine slots_reset

!!    run one batch: hru_control for every slot (in parallel with OpenMP), then the
!!    output and sums for every slot in command order
      subroutine hru_batch_run (ib)

      use hydrograph_module, only : ob, icmd, ht1, operator(*)
      use hru_module, only : ihru
      use time_module
      use basin_module, only : pco
      use carbon_module, only : org_con, org_allo

      integer, intent (in) :: ib      !none       |batch
      integer :: nslot                !none       |slots in the batch
      integer :: k                    !none       |slot
      integer :: iout                 !none       |counter
      integer :: ihtyp                !none       |hydrograph type out
      integer (kind=8) :: t0          !           |clock
      integer (kind=8) :: t1          !           |clock
      external :: hyddep_output, hydout_output

      if (stats_on) call system_clock (t0)
      nslot = batch_end(ib) - batch_beg(ib) + 1
      call slots_reset (nslot)
      call leftover_save_before

      !! org_con and org_allo contain carbon parameters as well as scratch
      !$omp parallel default (shared) private (k) copyin (org_con, org_allo)
      call dtbl_copies
      !$omp do schedule (dynamic, 1)
      do k = 1, nslot - 1
        call slot_run (ib, k)
      end do
      !$omp end do
      !$omp end parallel
      !! the last slot runs on this thread, so the module scratch left after the
      !! batch is the scratch the serial loop leaves after its last hru
      call slot_run (ib, nslot)
      islot = nslot + 1
      call leftover_load (islot)
      islot = 0
      if (stats_on) then
        call system_clock (t1)
        run_ticks = run_ticks + (t1 - t0)
        t0 = t1
      end if

      do k = 1, nslot
        icmd = batch_ob(batch_beg(ib) + k - 1)
        ihru = ob(icmd)%num
        call hru_sums (ihru)
        if (ob(icmd)%rcv_tot > 0) then
          ht1 = slot_ht1(k)
          call hyddep_output
        end if
        !print all outflow hydrographs
        if (time%yrs > pco%nyskip) then
          do iout = 1, ob(icmd)%src_tot
            ihtyp = ob(icmd)%ihtyp_out(iout)
            ht1 = ob(icmd)%frac_out(iout) * ob(icmd)%hd(ihtyp)
            call hydout_output (iout)
          end do
        end if
      end do
      if (stats_on) then
        call system_clock (t1)
        post_ticks = post_ticks + (t1 - t0)
      end if
      return
      end subroutine hru_batch_run

!!    run the hru in slot k of batch ib
      subroutine slot_run (ib, k)

      use hydrograph_module, only : ob, icmd, ht1
      use hru_module, only : ihru

      integer, intent (in) :: ib      !none       |batch
      integer, intent (in) :: k       !none       |slot
      integer :: iob                  !none       |last incoming object (not used)
      external :: obj_inflow, hru_control

      islot = k
      call leftover_mark
      icmd = batch_ob(batch_beg(ib) + k - 1)
      iob = 0
      call obj_inflow (iob)
      ihru = ob(icmd)%num
      call hru_control
      slot_ht1(k) = ht1
      call slot_finish (k)
      return
      end subroutine slot_run

!!    point d_tbl at land use decision table id
      subroutine dtbl_lum_point (id)

      use conditional_module, only : d_tbl, dtbl_lum

      integer, intent (in) :: id      !none       |table number in dtbl_lum

      if (islot == 0) then
        d_tbl => dtbl_lum(id)
      else if (lum_coupled(id)) then
        call ordered_gate ("dtbl_lum")
        d_tbl => dtbl_lum(id)
      else
        d_tbl => lum_copy(id)
      end if
      return
      end subroutine dtbl_lum_point

!!    point d_tbl at reservoir (wetland release) decision table id
      subroutine dtbl_res_point (id)

      use conditional_module, only : d_tbl, dtbl_res

      integer, intent (in) :: id      !none       |table number in dtbl_res

      if (islot == 0) then
        d_tbl => dtbl_res(id)
      else if (res_coupled(id)) then
        call ordered_gate ("dtbl_res")
        d_tbl => dtbl_res(id)
      else
        d_tbl => res_copy(id)
      end if
      return
      end subroutine dtbl_res_point

!!    refresh this thread's copies of the decision tables
      subroutine dtbl_copies

      use conditional_module, only : dtbl_lum, dtbl_res

      if (allocated (dtbl_lum)) lum_copy = dtbl_lum
      if (allocated (dtbl_res)) res_copy = dtbl_res
      return
      end subroutine dtbl_copies

!!    the leftover values of this thread
      subroutine leftover_get (r, n, arr)

      use hru_module, only : fertno3, fertnh3, fertorgn, fertsolp, fertorgp, enratio, etday, un2, up2, uapd

      real, dimension (nreal), intent (out) :: r          !   |leftover reals
      integer, dimension (nint), intent (out) :: n        !   |leftover integers
      real, dimension (:,:), intent (out) :: arr          !   |leftover arrays (array, plant)

      r = (/ fertno3, fertnh3, fertorgn, fertsolp, fertorgp, enratio, etday,           &
             lo_delg_p, lo_c_microb_sed, lo_gsi_adj, lo_yield /)
      n = (/ lo_j, lo_irrop, lo_idtill, lo_istr, lo_istr1 /)
      arr(1,:) = un2
      arr(2,:) = up2
      arr(3,:) = uapd
      return
      end subroutine leftover_get

!!    set the leftover values of this thread
      subroutine leftover_put (r, n, arr)

      use hru_module, only : fertno3, fertnh3, fertorgn, fertsolp, fertorgp, enratio, etday, un2, up2, uapd

      real, dimension (nreal), intent (in) :: r           !   |leftover reals
      integer, dimension (nint), intent (in) :: n         !   |leftover integers
      real, dimension (:,:), intent (in) :: arr           !   |leftover arrays (array, plant)

      fertno3 = r(1)
      fertnh3 = r(2)
      fertorgn = r(3)
      fertsolp = r(4)
      fertorgp = r(5)
      enratio = r(6)
      etday = r(7)
      lo_delg_p = r(8)
      lo_c_microb_sed = r(9)
      lo_gsi_adj = r(10)
      lo_yield = r(11)
      lo_j = n(1)
      lo_irrop = n(2)
      lo_idtill = n(3)
      lo_istr = n(4)
      lo_istr1 = n(5)
      un2 = arr(1,:)
      up2 = arr(2,:)
      uapd = arr(3,:)
      return
      end subroutine leftover_put

!!    keep the leftover values as it is before the batch (slot 0)
      subroutine leftover_save_before

      use mgt_operations_module, only : mgt
      use organic_mineral_mass_module, only : pl_mass_up

      slot_mgt(0) = mgt
      slot_mass_up(0) = pl_mass_up
      call leftover_get (slot_real(:,0), slot_int(:,0), slot_arr(:,:,0))
      slot_mgt_set(0) = .true.
      slot_mass_up_set(0) = .true.
      slot_real_set(:,0) = .true.
      slot_int_set(:,0) = .true.
      slot_arr_set(:,:,0) = .true.
      return
      end subroutine leftover_save_before

!!    mark the leftover values as not written at the start of a slot
      subroutine leftover_mark

      use mgt_operations_module, only : mgt
      use organic_mineral_mass_module, only : pl_mass_up

      real, dimension (nreal) :: r                        !   |leftover reals
      integer, dimension (nint) :: n                      !   |leftover integers
      real, dimension (narr,size(slot_arr,2)) :: arr      !   |leftover arrays

      mgt%husc = unset
      pl_mass_up%m = unset
      r = unset
      n = iunset
      arr = unset
      call leftover_put (r, n, arr)
      return
      end subroutine leftover_mark

!!    publish the leftover values slot k wrote
      subroutine leftover_save (k)

      use mgt_operations_module, only : mgt
      use organic_mineral_mass_module, only : pl_mass_up

      integer, intent (in) :: k                           !none   |slot
      real, dimension (nreal) :: r                        !   |leftover reals
      integer :: i                                        !none   |counter

      slot_mgt_set(k) = mgt%husc /= unset
      if (slot_mgt_set(k)) slot_mgt(k) = mgt
      slot_mass_up_set(k) = pl_mass_up%m /= unset
      if (slot_mass_up_set(k)) slot_mass_up(k) = pl_mass_up
      call leftover_get (r, slot_int(:,k), slot_arr(:,:,k))
      do i = 1, nreal
        !! values published at a milestone stay as they are (other slots may read them)
        if (any (mile_real == i .and. slot_mile(:,k))) cycle
        slot_real(i,k) = r(i)
        slot_real_set(i,k) = r(i) /= unset
      end do
      slot_int_set(:,k) = slot_int(:,k) /= iunset
      slot_arr_set(:,:,k) = slot_arr(:,:,k) /= unset
      return
      end subroutine leftover_save

!!    give the leftover values this slot has not written the value of the closest
!!    earlier slot that wrote it; every slot before k has finished
      subroutine leftover_load (k)

      use mgt_operations_module, only : mgt
      use organic_mineral_mass_module, only : pl_mass_up

      integer, intent (in) :: k                           !none   |slot
      real, dimension (nreal) :: r                        !   |leftover reals
      integer, dimension (nint) :: n                      !   |leftover integers
      real, dimension (narr,size(slot_arr,2)) :: arr      !   |leftover arrays
      integer :: i                                        !none   |counter
      integer :: ip                                       !none   |plant counter
      integer :: kk                                       !none   |earlier slot

      if (mgt%husc == unset) then
        do kk = k - 1, 0, -1
          if (slot_mgt_set(kk)) exit
        end do
        mgt = slot_mgt(kk)
      end if
      if (pl_mass_up%m == unset) then
        do kk = k - 1, 0, -1
          if (slot_mass_up_set(kk)) exit
        end do
        pl_mass_up = slot_mass_up(kk)
      end if
      call leftover_get (r, n, arr)
      do i = 1, nreal
        if (r(i) /= unset) cycle
        do kk = k - 1, 0, -1
          if (slot_real_set(i,kk)) exit
        end do
        r(i) = slot_real(i,kk)
      end do
      do i = 1, nint
        if (n(i) /= iunset) cycle
        do kk = k - 1, 0, -1
          if (slot_int_set(i,kk)) exit
        end do
        n(i) = slot_int(i,kk)
      end do
      do ip = 1, size(arr,2)
        do i = 1, narr
          if (arr(i,ip) /= unset) cycle
          do kk = k - 1, 0, -1
            if (slot_arr_set(i,ip,kk)) exit
          end do
          arr(i,ip) = slot_arr(i,ip,kk)
        end do
      end do
      call leftover_put (r, n, arr)
      return
      end subroutine leftover_load

!!    add the hru's part of the basin sums; called after hru_control in command order
      subroutine hru_sums (j)

      use basin_module, only : bsn_sedbud, pco
      use hru_module, only : sedyld
      use time_module

      integer, intent (in) :: j       !none       |hru number

      !! gwflow: add the hru's terms to the cell sums
      call gwdef_apply (j)
      if (wet_dep_on(j)) then
        bsn_sedbud%wet_dep_t = bsn_sedbud%wet_dep_t + wet_dep_d(j)
        wet_dep_on(j) = .false.
      end if
      !! don't sum during skip years
      if (time%yrs > pco%nyskip) then
        bsn_sedbud%upland_t = bsn_sedbud%upland_t + sedyld(j)
      end if
      return
      end subroutine hru_sums

!!    decision table classification
      logical function dtbl_coupled (dt)

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    returns .true. if evaluating the decision table with conditions and
!!    actions for one HRU can read or write state that another HRU in the
!!    batch can change, or writes to a file whose lines must stay in HRU
!!    order. Tables for which this is .true. have to run in serial order.
!!    Anything not recognised returns .true. (a false .true. costs speed;
!!    a false .false. breaks results).
!!    Classification: conditions.f90 and actions.f90 (58d6a1e).

!!    ~ ~ ~ NOTES ~ ~ ~
!!    - call after print.prt has been read: the 2612 rule uses pco%mgtout
!!    - reservoir tables (dtbl_res): the action type is "release" and the
!!      res_hydro option is in act(iac)%option; nonirr-h06, irr-h06 and
!!      hydrop write res_ob(jres) and return .true.
!!    - wetland_control adds to bsn_sedbud%wet_dep_t outside the table;
!!      HRUs with wetlands need a separate gate there
!!    - days_irr, irr_year, irr_demand and channel_flo are coupled only if
!!      an HRU is a water allocation POD (hru(j)%wallo_pod > 0, called at
!!      hru_control); the function cannot see that and returns .true.

      use basin_module, only : pco

      implicit none

      type (decision_table), intent (in) :: dt
      integer :: ic
      integer :: iac
      logical :: mgt_out
      logical :: ob_other
      character(len=25) :: cvar
      character(len=25) :: atyp
      character(len=25) :: aob
      character(len=40) :: aopt

      dtbl_coupled = .false.
      mgt_out = (pco%mgtout == "y")

      !! ======================================================================
      !! conditions (conditions.f90)
      !! ======================================================================
      if (dt%conds > 0) then
        if (.not. allocated(dt%cond)) then
          dtbl_coupled = .true.
          return
        end if
        if (size(dt%cond) < dt%conds) then
          dtbl_coupled = .true.
          return
        end if
      end if

      do ic = 1, dt%conds
        cvar = dt%cond(ic)%var
        ob_other = (dt%cond(ic)%ob_num /= 0)

        select case (trim(cvar))

        !! random numbers from the single stream rndseed_cond, and table state
        !! that persists across HRUs (days_prob, day_prev, prob_cum,
        !! hru_lu_cur, hru_ha_cur): conditions.f90:481-608
        case ("prob", "prob_unif", "prob_unif1", "prob_unif_lu")
          dtbl_coupled = .true.

        !! storage changed by irrigation withdrawals (actions.f90:267-318) and
        !! by water allocation (wallo_withdraw, wallo_return, wallo_pou_deliv):
        !! conditions.f90:741 (wtow_vol), 751-782 (vol)
        case ("wtow_vol", "vol")
          dtbl_coupled = .true.

        !! coupled only when an HRU is a water allocation POD
        !! (wallo_pou_deliv.f90:50,70,82 write pcom(j)%days_irr,
        !! hru(j)%irr_yr, irrig(j) of other HRUs); channel_flo reads the
        !! module scratch ht2, which wallo_control can change. The table
        !! cannot tell, so these return .true.
        case ("days_irr", "irr_year", "irr_demand", "channel_flo")
          dtbl_coupled = .true.

        !! time only, ob_num is not used: conditions.f90:263,313,319,343,349,355
        case ("day_start", "jday", "month", "year_cal", "year_seq", "year_start")
          continue

        !! data that no HRU changes during the batch, ob_num may be non-zero:
        !! ch_order sd_ch()%order, p_pet ob()%wst and wgn_pms(), slope,
        !! k_factor and hyd_soil_group static HRU and soil data (slope is
        !! changed only by calibration in serial code), aqu_dep dep_wt (only
        !! aqu_1d_control writes it), irr_pou_dmd pou()%demand (only
        !! wallo_start_day writes it), res_inflo ob()%hin of a reservoir (set
        !! by the reservoir command). aqu_dep and res_inflo also assume the
        !! batch keeps its place relative to the aquifer and reservoir commands.
        case ("ch_order", "p_pet", "slope", "k_factor", "hyd_soil_group",   &
              "aqu_dep", "irr_pou_dmd", "res_inflo")
          continue

        !! weather of the object's weather station, read only in the batch
        !! (climate_control and cli_precip_control run before it); iob is only
        !! set for ob = "hru" or "hlt" (conditions.f90:149-155), otherwise a
        !! stale SAVE value is used
        case ("phu_base0", "precip_cur", "precip_next")
          select case (trim(dt%cond(ic)%ob))
          case ("hru", "hlt")
            continue
          case default
            dtbl_coupled = .true.
          end select

        !! state of one HRU (pcom, pl_mass, soil, soil1, hru, hwb_d, wet_ob,
        !! wet): local when ob_num = 0, another HRU's state otherwise.
        !! tile_flo with ob_num /= 0 reads hwb_d() of another HRU, which is
        !! written at the end of its hru_control (hru_control.f90:765)
        case ("w_stress", "n_stress", "p_stress", "phu_plant", "plant_gro",  &
              "plant_name_gro", "days_plant", "days_harv", "days_act",       &
              "soil_water", "year_rot", "year_gro", "cur_yrs_mat",           &
              "biomass", "leaf_area", "ground_cov", "p_factor",              &
              "p_lab_150", "soil_temp2", "soil_carbon", "tile_drained",      &
              "tile_flo", "land_use", "cal_group", "tillage", "plant_com",   &
              "wet_depth", "weirh", "vol_wet")
          if (ob_other) dtbl_coupled = .true.

        !! unknown condition variable
        case default
          dtbl_coupled = .true.

        end select

        if (dtbl_coupled) return
      end do

      !! ======================================================================
      !! actions (actions.f90)
      !! ======================================================================
      if (dt%acts > 0) then
        if (.not. allocated(dt%act)) then
          dtbl_coupled = .true.
          return
        end if
        if (size(dt%act) < dt%acts) then
          dtbl_coupled = .true.
          return
        end if
      end if

      do iac = 1, dt%acts
        atyp = dt%act(iac)%typ
        aopt = dt%act(iac)%option
        aob = dt%act(iac)%ob
        ob_other = (dt%act(iac)%ob_num /= 0)

        select case (trim(atyp))

        !! writes other objects or shared and basin state: mallo, res_ob, pou,
        !! module dmd_m3, trn_fr, basin crop yield and plcal sums,
        !! hru(hru_rcv)%sb%inflo, ru_elem, lu_chg_mx, snow_chg_mx,
        !! basin_plants, plts_bsn, sd_ch. grow_init and grow_end index hlt(),
        !! the hru_lte array, with an HRU number. irr_wallo has no case in
        !! actions.f90 but belongs to water allocation.
        !! actions.f90:112-122, 171-237, 495-577, 611-694, 811-865, 979-1045,
        !! 1218-1222, 770-808
        case ("manure_demand", "res_irr_dmd", "duty_reset", "wallo_dmd",      &
              "duty_fr", "tileflo_control",                                   &
              "hru_fr_update", "lu_change", "snow_change", "chan_change",     &
              "grow_init", "grow_end", "irr_wallo")
          dtbl_coupled = .true.

        !! 3612 is written unconditionally, so file lines would interleave:
        !! actions.f90:1060,1075,1092,1109,1130,1149,1168,1188,1196,1214
        case ("p_factor", "contour", "stripcrop", "terrace", "tile_install",  &
              "septic_install", "fstrip_install", "grassww_install",          &
              "grassww_uninstall", "user_def_bmp")
          dtbl_coupled = .true.

        !! act%ob_num is the source object (aqu, cha, sdc, res), not a target
        !! HRU. actions calls the gate before it takes water from a source, and
        !! harvest and harvest_kill call it before they add to the basin yields,
        !! so the rest of the table can run on the thread copy
        case ("irrigate")
          if (mgt_out) dtbl_coupled = .true.
        case ("harvest", "harvest_kill")
          if (mgt_out .or. ob_other) dtbl_coupled = .true.

        !! act%ob_num is the fertilizer slot, not an object: actions.f90:405
        case ("fert_future")
          continue

        !! weir height: option "wet" writes wet_ob(j) only; any other option
        !! writes res_ob(j)%weir_hgt, which is reservoir j (actions.f90:917-934)
        case ("weir_height")
          if (trim(aopt) /= "wet") dtbl_coupled = .true.
          if (ob_other) dtbl_coupled = .true.
          if (mgt_out) dtbl_coupled = .true.

        !! writes only the HRU j = act%ob_num (0 = current HRU) and unit 2612
        !! when pco%mgtout == "y"; another HRU when ob_num /= 0.
        !! actions.f90:126-168, 344-399, 408-429, 432-492, 580-608, 732-752,
        !! 870-884, 887-914, 937-976, 1225-1244, 1247-1264
        case ("irr_demand", "fertilize", "manure", "till", "tillage", "plant",   &
              "kill", "pest_apply", "tiledep_control", "impound_off",            &
              "impound_on", "puddle", "burn", "cn_update")
          if (ob_other) dtbl_coupled = .true.
          if (mgt_out) dtbl_coupled = .true.
          !! transplanting reads pl_mass_up before writing it (leftover values)
          if (trim(atyp) == "plant" .and. allocated(dt%act_app)) then
            if (dt%act_app(iac) > 0) dtbl_coupled = .true.
          end if

        !! writes only the HRU j = act%ob_num (0 = current HRU), no file
        !! actions.f90:697-729, 754-767, 1266-1283
        case ("rot_reset", "harv_reset", "kill_reset", "plant_reset",        &
              "irr_reset", "pheno_reset")
          if (ob_other) dtbl_coupled = .true.

        !! graze sets the module variable graze, which every grazing hru reads
        case ("graze")
          dtbl_coupled = .true.

        !! empty case: actions.f90:1286
        case ("herd")
          continue

        !! reservoir and wetland release tables (res_hydro.f90:103-343): the
        !! options below touch ht1, ht2, wbody and wet_ob(jres) only; nonirr-h06,
        !! irr-h06 and hydrop write res_ob(jres) (reservoir number = HRU number
        !! when called from wetland_control), and so does any unknown option
        case ("release")
          select case (trim(aopt))
          case ("rate", "rate_pct", "inflo_rate", "inflo_frac", "ab_emer",    &
                "days", "dyrt", "dyrt1", "inflo_targ", "weir", "meas",        &
                "natlake")
            continue
          case default
            dtbl_coupled = .true.
          end select

        !! unknown action type
        case default
          dtbl_coupled = .true.

        end select

        if (dtbl_coupled) return
      end do

      return
      end function dtbl_coupled

      end module hru_batch_module
