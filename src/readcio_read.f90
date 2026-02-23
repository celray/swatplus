      subroutine readcio_read

       use input_file_module
       use output_path_module
       use utils, only: to_lower

       implicit none

       character (len=25) :: section_name = ""
       character (len=256) :: out_path_value = ""
       character (len=1024) :: line_buffer = ""
       character (len=1024) :: line_rest = ""
       character (len=256) :: tokens(50)
       integer :: ntokens = 0
       integer :: eof = 0
       integer :: idx = 0
       logical :: i_exist
       logical :: found_io_path = .false.

       eof = 0

       !! read file.cio
       inquire (file="file.cio", exist=i_exist)
       if (.not. i_exist) then
         write (*,*) 'ERROR: file.cio not found'
         stop 1
       end if

       open (107, file="file.cio")
       read (107, '(A)', iostat=eof) line_buffer    !! skip header line

       !! Section-name-based reader for the redesigned file.cio format
       do
         read (107, '(A)', iostat=eof) line_buffer
         if (eof < 0) exit

         !! skip blank lines
         line_buffer = adjustl(line_buffer)
         if (len_trim(line_buffer) == 0) cycle

         !! extract section name (first whitespace-delimited token)
         idx = index(trim(line_buffer), ' ')
         if (idx == 0) then
           section_name = trim(line_buffer)
           line_rest = ""
         else
           section_name = line_buffer(1:idx-1)
           line_rest = adjustl(line_buffer(idx+1:))
         end if

         call parse_line_tokens(line_rest, tokens, ntokens)

         section_name = to_lower(trim(section_name))

         select case (trim(section_name))

         !! simulation: time.sim  print.prt  object.prt  object.cnt
         case ('simulation')
           if (ntokens >= 1) in_sim%time       = tokens(1)
           if (ntokens >= 2) in_sim%prt        = tokens(2)
           if (ntokens >= 3) in_sim%object_prt = tokens(3)
           if (ntokens >= 4) in_sim%object_cnt = tokens(4)

         !! basin: codes.bsn  parameters.bsn
         case ('basin')
           if (ntokens >= 1) in_basin%codes_bas = tokens(1)
           if (ntokens >= 2) in_basin%parms_bas = tokens(2)

         !! climate: weather_sta.cli  weather_wgn.cli  pet.cli  pcp.cli  tmp.cli  slr.cli  hmd.cli  wnd.cli  atmodep.cli  atmosalt.cli  atmocs.cli
         case ('climate')
           if (ntokens >= 1)  in_cli%weat_sta  = tokens(1)
           if (ntokens >= 2)  in_cli%weat_wgn  = tokens(2)
           if (ntokens >= 3)  in_cli%pet_cli   = tokens(3)
           if (ntokens >= 4)  in_cli%pcp_cli   = tokens(4)
           if (ntokens >= 5)  in_cli%tmp_cli   = tokens(5)
           if (ntokens >= 6)  in_cli%slr_cli   = tokens(6)
           if (ntokens >= 7)  in_cli%hmd_cli   = tokens(7)
           if (ntokens >= 8)  in_cli%wnd_cli   = tokens(8)
           if (ntokens >= 9)  in_cli%atmo_cli  = tokens(9)
           if (ntokens >= 10) in_cli%atmo_salt = tokens(10)
           if (ntokens >= 11) in_cli%atmo_cs   = tokens(11)

         !! connect: hru.con  hru_lte.con  rout_unit.con  gwflow.con  aquifer.con  aquifer2d.con  channel.con  reservoir.con  recall.con  exco.con  delratio.con  outlet.con  chandeg.con
         case ('connect')
           if (ntokens >= 1)  in_con%hru_con     = tokens(1)
           if (ntokens >= 2)  in_con%hruez_con   = tokens(2)
           if (ntokens >= 3)  in_con%ru_con      = tokens(3)
           if (ntokens >= 4)  in_con%gwflow_con  = tokens(4)
           if (ntokens >= 5)  in_con%aqu_con     = tokens(5)
           if (ntokens >= 6)  in_con%aqu2d_con   = tokens(6)
           if (ntokens >= 7)  in_con%chan_con     = tokens(7)
           if (ntokens >= 8)  in_con%res_con     = tokens(8)
           if (ntokens >= 9)  in_con%rec_con     = tokens(9)
           if (ntokens >= 10) in_con%exco_con    = tokens(10)
           if (ntokens >= 11) in_con%delr_con    = tokens(11)
           if (ntokens >= 12) in_con%out_con     = tokens(12)
           if (ntokens >= 13) in_con%chandeg_con = tokens(13)

         !! channel: initial.cha  channel.cha  hydrology.cha  sediment.cha  nutrients.cha  channel_lte.cha  hyd_sed_lte.cha  temperature.cha  sed_nut.cha  element.ccu
         case ('channel')
           if (ntokens >= 1)  in_cha%init        = tokens(1)
           if (ntokens >= 2)  in_cha%dat         = tokens(2)
           if (ntokens >= 3)  in_cha%hyd         = tokens(3)
           if (ntokens >= 4)  in_cha%sed         = tokens(4)
           if (ntokens >= 5)  in_cha%nut         = tokens(5)
           if (ntokens >= 6)  in_cha%chan_ez     = tokens(6)
           if (ntokens >= 7)  in_cha%hyd_sed     = tokens(7)
           if (ntokens >= 8)  in_cha%temp        = tokens(8)
           if (ntokens >= 9)  in_cha%sed_nut     = tokens(9)
           if (ntokens >= 10) in_cha%element_ccu = tokens(10)

         !! reservoir: initial.res  reservoir.res  hydrology.res  sediment.res  nutrients.res  weir.res  wetland.wet  hydrology.wet
         case ('reservoir')
           if (ntokens >= 1) in_res%init_res = tokens(1)
           if (ntokens >= 2) in_res%res      = tokens(2)
           if (ntokens >= 3) in_res%hyd_res  = tokens(3)
           if (ntokens >= 4) in_res%sed_res  = tokens(4)
           if (ntokens >= 5) in_res%nut_res  = tokens(5)
           if (ntokens >= 6) in_res%weir_res = tokens(6)
           if (ntokens >= 7) in_res%wet       = tokens(7)
           if (ntokens >= 8) in_res%hyd_wet   = tokens(8)
           if (ntokens >= 9) in_res%res_conds = tokens(9)

         !! routing_unit: rout_unit.def  rout_unit.ele  rout_unit.rtu  rout_unit.dr
         case ('routing_unit')
           if (ntokens >= 1) in_ru%ru_def = tokens(1)
           if (ntokens >= 2) in_ru%ru_ele = tokens(2)
           if (ntokens >= 3) in_ru%ru     = tokens(3)
           if (ntokens >= 4) in_ru%ru_dr  = tokens(4)

         !! hru: hru_data.hru  hru_lte.hru
         case ('hru')
           if (ntokens >= 1) in_hru%hru_data = tokens(1)
           if (ntokens >= 2) in_hru%hru_ez   = tokens(2)

         !! exco: exco.exc  exco_om.exc  exco_pest.exc  exco_path.exc  exco_hmet.exc  exco_salt.exc
         case ('exco')
           if (ntokens >= 1) in_exco%exco = tokens(1)
           if (ntokens >= 2) in_exco%om   = tokens(2)
           if (ntokens >= 3) in_exco%pest = tokens(3)
           if (ntokens >= 4) in_exco%path = tokens(4)
           if (ntokens >= 5) in_exco%hmet = tokens(5)
           if (ntokens >= 6) in_exco%salt = tokens(6)

         !! recall: recall.rec  recall.slt  recall.cs
         case ('recall')
           if (ntokens >= 1) in_rec%recall_rec = tokens(1)
           if (ntokens >= 2) in_rec%recall_slt = tokens(2)
           if (ntokens >= 3) in_rec%recall_cs  = tokens(3)
           if (ntokens >= 4) in_rec%recall_db  = tokens(4)
           if (ntokens >= 5) in_rec%pest_com   = tokens(5)

         !! dr: delratio.del  dr_om.del  dr_pest.del  dr_path.del  dr_hmet.del  dr_salt.del
         case ('dr')
           if (ntokens >= 1) in_delr%del_ratio = tokens(1)
           if (ntokens >= 2) in_delr%om        = tokens(2)
           if (ntokens >= 3) in_delr%pest      = tokens(3)
           if (ntokens >= 4) in_delr%path      = tokens(4)
           if (ntokens >= 5) in_delr%hmet      = tokens(5)
           if (ntokens >= 6) in_delr%salt      = tokens(6)

         !! aquifer: initial.aqu  aquifer.aqu  gwflow.aqu
         case ('aquifer')
           if (ntokens >= 1) in_aqu%init   = tokens(1)
           if (ntokens >= 2) in_aqu%aqu    = tokens(2)
           if (ntokens >= 3) in_aqu%gwflow = tokens(3)

         !! herd: animal.hrd  herd.hrd  ranch.hrd
         case ('herd')
           if (ntokens >= 1) in_herd%animal = tokens(1)
           if (ntokens >= 2) in_herd%herd   = tokens(2)
           if (ntokens >= 3) in_herd%ranch  = tokens(3)

         !! link: chan_surf.lin  aqu_cha.lin
         case ('link')
           if (ntokens >= 1) in_link%chan_surf = tokens(1)
           if (ntokens >= 2) in_link%aqu_cha  = tokens(2)

         !! hydrology: hydrology.hyd  topography.hyd  field.fld
         case ('hydrology')
           if (ntokens >= 1) in_hyd%hydrol_hyd = tokens(1)
           if (ntokens >= 2) in_hyd%topogr_hyd = tokens(2)
           if (ntokens >= 3) in_hyd%field_fld  = tokens(3)

         !! structural: tiledrain.str  septic.str  filterstrip.str  grassedww.str  bmpuser.str  satbuffer.str
         case ('structural')
           if (ntokens >= 1) in_str%tiledrain_str  = tokens(1)
           if (ntokens >= 2) in_str%septic_str     = tokens(2)
           if (ntokens >= 3) in_str%fstrip_str     = tokens(3)
           if (ntokens >= 4) in_str%grassww_str    = tokens(4)
           if (ntokens >= 5) in_str%bmpuser_str    = tokens(5)
           if (ntokens >= 6) in_str%satbuffer_str  = tokens(6)

         !! hru_parm_db: plants.plt  fertilizer.frt  tillage.til  pesticide.pes  metabolite.pes  pathogens.pth  metals.mtl  salts.slt  urban.urb  septic.sep  snow.sno
         case ('hru_parm_db')
           if (ntokens >= 1)  in_parmdb%plants_plt = tokens(1)
           if (ntokens >= 2)  in_parmdb%fert_frt   = tokens(2)
           if (ntokens >= 3)  in_parmdb%till_til   = tokens(3)
           if (ntokens >= 4)  in_parmdb%pest       = tokens(4)
           if (ntokens >= 5)  in_parmdb%metabolite = tokens(5)
           if (ntokens >= 6)  in_parmdb%pathcom_db = tokens(6)
           if (ntokens >= 7)  in_parmdb%hmetcom_db = tokens(7)
           if (ntokens >= 8)  in_parmdb%saltcom_db = tokens(8)
           if (ntokens >= 9)  in_parmdb%urban_urb  = tokens(9)
           if (ntokens >= 10) in_parmdb%septic_sep = tokens(10)
           if (ntokens >= 11) in_parmdb%snow       = tokens(11)

         !! ops: harv.ops  graze.ops  irr.ops  chem_app.ops  fire.ops  sweep.ops  puddle.ops  transplant.ops
         case ('ops')
           if (ntokens >= 1) in_ops%harv_ops       = tokens(1)
           if (ntokens >= 2) in_ops%graze_ops      = tokens(2)
           if (ntokens >= 3) in_ops%irr_ops        = tokens(3)
           if (ntokens >= 4) in_ops%chem_ops       = tokens(4)
           if (ntokens >= 5) in_ops%fire_ops       = tokens(5)
           if (ntokens >= 6) in_ops%sweep_ops      = tokens(6)
           if (ntokens >= 7) in_ops%puddle_ops     = tokens(7)
           if (ntokens >= 8) in_ops%transplant_ops = tokens(8)

         !! lum: landuse.lum  management.sch  cntable.lum  cons_practice.lum  ovn_table.lum
         case ('lum')
           if (ntokens >= 1) in_lum%landuse_lum    = tokens(1)
           if (ntokens >= 2) in_lum%management_sch = tokens(2)
           if (ntokens >= 3) in_lum%cntable_lum    = tokens(3)
           if (ntokens >= 4) in_lum%cons_prac_lum  = tokens(4)
           if (ntokens >= 5) in_lum%ovn_lum        = tokens(5)

         !! calibration: cal_parms.cal  calibration.cal  codes.sft  wb_parms.sft  water_balance.sft  ch_sed_budget.sft  ch_sed_parms.sft  plant_parms.sft  plant_gro.sft
         case ('calibration')
           if (ntokens >= 1) in_chg%cal_parms          = tokens(1)
           if (ntokens >= 2) in_chg%cal_upd            = tokens(2)
           if (ntokens >= 3) in_chg%codes_sft          = tokens(3)
           if (ntokens >= 4) in_chg%wb_parms_sft       = tokens(4)
           if (ntokens >= 5) in_chg%water_balance_sft  = tokens(5)
           if (ntokens >= 6) in_chg%ch_sed_budget_sft  = tokens(6)
           if (ntokens >= 7) in_chg%ch_sed_parms_sft   = tokens(7)
           if (ntokens >= 8) in_chg%plant_parms_sft    = tokens(8)
           if (ntokens >= 9) in_chg%plant_gro_sft      = tokens(9)

         !! init: plant.ini  soil_plant.ini  om_water.ini  pest_hru.ini  pest_water.ini  path_hru.ini  path_water.ini  hmet_hru.ini  hmet_water.ini  salt_hru.ini  salt_water.ini
         case ('init')
           if (ntokens >= 1)  in_init%plant          = tokens(1)
           if (ntokens >= 2)  in_init%soil_plant_ini = tokens(2)
           if (ntokens >= 3)  in_init%om_water       = tokens(3)
           if (ntokens >= 4)  in_init%pest_soil      = tokens(4)
           if (ntokens >= 5)  in_init%pest_water     = tokens(5)
           if (ntokens >= 6)  in_init%path_soil      = tokens(6)
           if (ntokens >= 7)  in_init%path_water     = tokens(7)
           if (ntokens >= 8)  in_init%hmet_soil      = tokens(8)
           if (ntokens >= 9)  in_init%hmet_water     = tokens(9)
           if (ntokens >= 10) in_init%salt_soil      = tokens(10)
           if (ntokens >= 11) in_init%salt_water     = tokens(11)

         !! soils: soils.sol  nutrients.sol  soils_lte.sol
         case ('soils')
           if (ntokens >= 1) in_sol%soils_sol   = tokens(1)
           if (ntokens >= 2) in_sol%nut_sol     = tokens(2)
           if (ntokens >= 3) in_sol%lte_sol     = tokens(3)
           if (ntokens >= 4) in_sol%lyr_depths  = tokens(4)

         !! decision_table: lum.dtl  res_rel.dtl  scen_lu.dtl  flo_con.dtl
         case ('decision_table')
           if (ntokens >= 1) in_cond%dtbl_lum  = tokens(1)
           if (ntokens >= 2) in_cond%dtbl_res  = tokens(2)
           if (ntokens >= 3) in_cond%dtbl_scen = tokens(3)
           if (ntokens >= 4) in_cond%dtbl_flo  = tokens(4)

         !! regions: ls_unit.ele  ls_unit.def  ls_reg.ele  ls_reg.def  ls_cal.reg  ch_catunit.ele  ch_catunit.def  ch_reg.def  aqu_catunit.ele  aqu_catunit.def  aqu_reg.def  res_catunit.ele  res_catunit.def  res_reg.def  rec_catunit.ele  rec_catunit.def  rec_reg.def
         case ('regions')
           if (ntokens >= 1)  in_regs%ele_lsu     = tokens(1)
           if (ntokens >= 2)  in_regs%def_lsu     = tokens(2)
           if (ntokens >= 3)  in_regs%ele_reg     = tokens(3)
           if (ntokens >= 4)  in_regs%def_reg     = tokens(4)
           if (ntokens >= 5)  in_regs%cal_lcu     = tokens(5)
           if (ntokens >= 6)  in_regs%ele_cha     = tokens(6)
           if (ntokens >= 7)  in_regs%def_cha     = tokens(7)
           if (ntokens >= 8)  in_regs%def_cha_reg = tokens(8)
           if (ntokens >= 9)  in_regs%ele_aqu     = tokens(9)
           if (ntokens >= 10) in_regs%def_aqu     = tokens(10)
           if (ntokens >= 11) in_regs%def_aqu_reg = tokens(11)
           if (ntokens >= 12) in_regs%ele_res     = tokens(12)
           if (ntokens >= 13) in_regs%def_res     = tokens(13)
           if (ntokens >= 14) in_regs%def_res_reg = tokens(14)
           if (ntokens >= 15) in_regs%ele_psc     = tokens(15)
           if (ntokens >= 16) in_regs%def_psc     = tokens(16)
           if (ntokens >= 17) in_regs%def_psc_reg = tokens(17)

         !! carbon: basins.cbn  coefficients.cbn  co2_yr.cbn
         case ('carbon')
           if (ntokens >= 1) in_carbon%basins_cbn       = tokens(1)
           if (ntokens >= 2) in_carbon%coefficients_cbn = tokens(2)
           if (ntokens >= 3) in_carbon%co2_yr_cbn       = tokens(3)

         !! salt: initial.slt  channel.slt  hru.slt  fertilizer.slt  irrigation.slt  plants.slt  road.slt  urban.slt  uptake.slt  reservoir.slt
         case ('salt')
           if (ntokens >= 1)  in_salt%init_slt    = tokens(1)
           if (ntokens >= 2)  in_salt%channel_slt = tokens(2)
           if (ntokens >= 3)  in_salt%hru_slt     = tokens(3)
           if (ntokens >= 4)  in_salt%fert_slt    = tokens(4)
           if (ntokens >= 5)  in_salt%irr_slt     = tokens(5)
           if (ntokens >= 6)  in_salt%plants_slt  = tokens(6)
           if (ntokens >= 7)  in_salt%road_slt    = tokens(7)
           if (ntokens >= 8)  in_salt%urban_slt   = tokens(8)
           if (ntokens >= 9)  in_salt%uptake_slt  = tokens(9)
           if (ntokens >= 10) in_salt%res_slt     = tokens(10)

         !! constituents: constituents.cs  initial.cs  channel.cs  hru.cs  fertilizer.cs  irrigation.cs  plants_boron.cs  reactions.cs  uptake.cs  urban.cs  reservoir.cs  streamobs.cs  initial_aqu.cs  initial_cha.cs  reservoir_cs.cs  wetland.cs  nutrients.rte
         case ('constituents')
           if (ntokens >= 1)  in_constit%cs_db           = tokens(1)
           if (ntokens >= 1)  in_sim%cs_db               = tokens(1)   !! also set in_sim%cs_db for existing code that references it
           if (ntokens >= 2)  in_constit%init_cs         = tokens(2)
           if (ntokens >= 3)  in_constit%channel_cs      = tokens(3)
           if (ntokens >= 4)  in_constit%hru_cs          = tokens(4)
           if (ntokens >= 5)  in_constit%fert_cs         = tokens(5)
           if (ntokens >= 6)  in_constit%irr_cs          = tokens(6)
           if (ntokens >= 7)  in_constit%plants_boron_cs = tokens(7)
           if (ntokens >= 8)  in_constit%reactions_cs    = tokens(8)
           if (ntokens >= 9)  in_constit%uptake_cs       = tokens(9)
           if (ntokens >= 10) in_constit%urban_cs        = tokens(10)
           if (ntokens >= 11) in_constit%res_cs          = tokens(11)
           if (ntokens >= 12) in_constit%streamobs_cs    = tokens(12)
           if (ntokens >= 13) in_constit%init_aqu_cs     = tokens(13)
           if (ntokens >= 14) in_constit%init_cha_cs     = tokens(14)
           if (ntokens >= 15) in_constit%reservoir_cs    = tokens(15)
           if (ntokens >= 16) in_constit%wetland_cs      = tokens(16)
           if (ntokens >= 17) in_constit%nut_rte         = tokens(17)

         !! manure: manure.frt  manure_allo.mnu
         case ('manure')
           if (ntokens >= 1) in_manure%manure_frt  = tokens(1)
           if (ntokens >= 2) in_manure%manure_allo = tokens(2)

         !! water_allocation: water_allocation.wro  water_pipe.wal  water_tower.wal  water_use.wal  water_treat.wal  om_treat.wal  om_use.wal
         case ('water_allocation')
           if (ntokens >= 1) in_wal%wro         = tokens(1)
           if (ntokens >= 2) in_wal%water_pipe  = tokens(2)
           if (ntokens >= 3) in_wal%water_tower = tokens(3)
           if (ntokens >= 4) in_wal%water_use   = tokens(4)
           if (ntokens >= 5) in_wal%water_treat = tokens(5)
           if (ntokens >= 6) in_wal%om_treat     = tokens(6)
           if (ntokens >= 7) in_wal%om_use       = tokens(7)
           if (ntokens >= 8) in_wal%outside_src  = tokens(8)
           if (ntokens >= 9) in_wal%om_osrc      = tokens(9)

         !! update: scen_dtl.upd
         case ('update')
           if (ntokens >= 1) in_upd%scen_dtl = tokens(1)

         !! io_path: pcp_path  tmp_path  slr_path  hmd_path  wnd_path  pet_path  out_path
         case ('io_path')
           found_io_path = .true.
           if (ntokens >= 1) in_path_pcp%pcp  = tokens(1)
           if (ntokens >= 2) in_path_tmp%tmp  = tokens(2)
           if (ntokens >= 3) in_path_slr%slr  = tokens(3)
           if (ntokens >= 4) in_path_hmd%hmd  = tokens(4)
           if (ntokens >= 5) in_path_wnd%wnd  = tokens(5)
           if (ntokens >= 6) in_path_pet%peti = tokens(6)
           if (ntokens >= 7) out_path_value   = tokens(7)

         case default
           continue

         end select
       end do

       close (107)

       !! Check for old-format file.cio
       if (.not. found_io_path) then
         write (*,*) ''
         write (*,*) '! Error - file.cio is in an old format that is no longer supported.'
         write (*,*) '!         Please run the file.cio conversion tool from the SWAT+ website'
         write (*,*) '!         in the TxtInOut directory to update it to the current format.'
         stop 1
       end if

       !! Initialize output path (will use current dir if null/empty)
       call init_output_path(out_path_value)

       return

       contains

       !! parse whitespace-delimited tokens from a line
       subroutine parse_line_tokens(line, toks, count)
         character(len=*), intent(in) :: line
         character(len=256), intent(out) :: toks(50)
         integer, intent(out) :: count

         character(len=1024) :: buf
         integer :: i, j, slen

         count = 0
         toks = ""
         buf = adjustl(line)
         slen = len_trim(buf)
         if (slen == 0) return

         i = 1
         do while (i <= slen .and. count < 50)
           !! skip whitespace
           do while (i <= slen .and. (buf(i:i) == ' ' .or. buf(i:i) == char(9)))
             i = i + 1
           end do
           if (i > slen) exit

           !! found start of token
           j = i
           do while (j <= slen .and. buf(j:j) /= ' ' .and. buf(j:j) /= char(9))
             j = j + 1
           end do

           count = count + 1
           toks(count) = buf(i:j-1)
           i = j
         end do
       end subroutine parse_line_tokens

      end subroutine readcio_read
