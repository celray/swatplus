      module input_file_module

      implicit none

!! simulation
      type input_sim
        character(len=25) :: time = ""
        character(len=25) :: prt = ""
        character(len=25) :: object_prt = ""
        character(len=25) :: object_cnt = ""
        character(len=25) :: cs_db = ""           !! kept for code that still references in_sim%cs_db; canonical location is in_constit%cs_db
      end type input_sim
      type (input_sim) :: in_sim

!! basin
      type input_basin
       character(len=25) :: codes_bas = ""
       character(len=25) :: parms_bas = ""
      end type input_basin
      type (input_basin) :: in_basin

!! climate
      type input_cli
       character(len=25) :: weat_sta = ""
       character(len=25) :: weat_wgn = ""
       character(len=25) :: pet_cli = ""
       character(len=25) :: pcp_cli = ""
       character(len=25) :: tmp_cli = ""
       character(len=25) :: slr_cli = ""
       character(len=25) :: hmd_cli = ""
       character(len=25) :: wnd_cli = ""
       character(len=25) :: atmo_cli = ""
       character(len=25) :: atmo_salt = ""
       character(len=25) :: atmo_cs = ""
      end type input_cli
      type (input_cli) :: in_cli

!! connect
      type input_con
       character(len=25) :: hru_con = ""
       character(len=25) :: hruez_con = ""
       character(len=25) :: ru_con = ""
       character(len=25) :: gwflow_con = ""
       character(len=25) :: aqu_con = ""
       character(len=25) :: aqu2d_con = ""
       character(len=25) :: chan_con = ""
       character(len=25) :: res_con = ""
       character(len=25) :: rec_con = ""
       character(len=25) :: exco_con = ""
       character(len=25) :: delr_con = ""
       character(len=25) :: out_con = ""
       character(len=25) :: chandeg_con = ""
      end type input_con
      type (input_con) :: in_con

!! channel
      type input_cha
       character(len=25) :: init = ""
       character(len=25) :: dat = ""
       character(len=25) :: hyd = ""
       character(len=25) :: sed = ""
       character(len=25) :: nut = ""
       character(len=25) :: chan_ez = ""
       character(len=25) :: hyd_sed = ""
       character(len=25) :: temp = ""
       character(len=25) :: sed_nut = ""
       character(len=25) :: element_ccu = ""
      end type input_cha
      type (input_cha) :: in_cha

!! reservoir
      type input_res
       character(len=25) :: init_res = ""
       character(len=25) :: res = ""
       character(len=25) :: hyd_res = ""
       character(len=25) :: sed_res = ""
       character(len=25) :: nut_res = ""
       character(len=25) :: weir_res = ""
       character(len=25) :: wet = ""
       character(len=25) :: hyd_wet = ""
       character(len=25) :: res_conds = ""
      end type input_res
      type (input_res) :: in_res

!! routing unit
      type input_ru
       character(len=25) :: ru_def = ""
       character(len=25) :: ru_ele = ""
       character(len=25) :: ru = ""
       character(len=25) :: ru_dr = ""
      end type input_ru
      type (input_ru) :: in_ru

!! HRU
      type input_hru
       character(len=25) :: hru_data = ""
       character(len=25) :: hru_ez = ""
      end type input_hru
      type (input_hru) :: in_hru

!! exco
      type input_exco
       character(len=25) :: exco = ""
       character(len=25) :: om = ""
       character(len=25) :: pest = ""
       character(len=25) :: path = ""
       character(len=25) :: hmet = ""
       character(len=25) :: salt = ""
      end type input_exco
      type (input_exco) :: in_exco

!! recall
      type input_rec
       character(len=25) :: recall_rec = ""
       character(len=25) :: recall_slt = ""
       character(len=25) :: recall_cs = ""
       character(len=25) :: recall_db = ""
       character(len=25) :: pest_com = ""
      end type input_rec
      type (input_rec) :: in_rec

!! delivery ratio
      type input_delr
       character(len=25) :: del_ratio = ""
       character(len=25) :: om = ""
       character(len=25) :: pest = ""
       character(len=25) :: path = ""
       character(len=25) :: hmet = ""
       character(len=25) :: salt = ""
      end type input_delr
      type (input_delr) :: in_delr

!! aquifer
      type input_aqu
       character(len=25) :: init = ""
       character(len=25) :: aqu = ""
       character(len=25) :: gwflow = ""
      end type input_aqu
      type (input_aqu) :: in_aqu

!! herd
      type input_herd
        character(len=25) :: animal = ""
        character(len=25) :: herd = ""
        character(len=25) :: ranch = ""
      end type input_herd
      type (input_herd) :: in_herd

!! link
      type input_link
       character(len=25) :: chan_surf = ""
       character(len=25) :: aqu_cha = ""
      end type input_link
      type (input_link) :: in_link

!! hydrology
      type input_hydrology
       character(len=25) :: hydrol_hyd = ""
       character(len=25) :: topogr_hyd = ""
       character(len=25) :: field_fld = ""
      end type input_hydrology
      type (input_hydrology) :: in_hyd

!! structural
      type input_structural
       character(len=25) :: tiledrain_str = ""
       character(len=25) :: septic_str = ""
       character(len=25) :: fstrip_str = ""
       character(len=25) :: grassww_str = ""
       character(len=25) :: bmpuser_str = ""
       character(len=25) :: satbuffer_str = ""
      end type input_structural
      type (input_structural) :: in_str

!! HRU parameter databases
      type input_parameter_databases
       character(len=25) :: plants_plt = ""
       character(len=25) :: fert_frt = ""
       character(len=25) :: till_til = ""
       character(len=25) :: pest = ""
       character(len=25) :: metabolite = ""
       character(len=25) :: pathcom_db = ""
       character(len=25) :: hmetcom_db = ""
       character(len=25) :: saltcom_db = ""
       character(len=25) :: urban_urb = ""
       character(len=25) :: septic_sep = ""
       character(len=25) :: snow = ""
      end type input_parameter_databases
      type (input_parameter_databases) :: in_parmdb

!! operation scheduling
      type input_ops
       character(len=25) :: harv_ops = ""
       character(len=25) :: graze_ops = ""
       character(len=25) :: irr_ops = ""
       character(len=25) :: chem_ops = ""
       character(len=25) :: fire_ops = ""
       character(len=25) :: sweep_ops = ""
       character(len=25) :: puddle_ops = ""
       character(len=25) :: transplant_ops = ""
      end type input_ops
      type (input_ops) :: in_ops

!! land use management
      type input_lum
       character(len=25) :: landuse_lum = ""
       character(len=25) :: management_sch = ""
       character(len=25) :: cntable_lum = ""
       character(len=25) :: cons_prac_lum = ""
       character(len=25) :: ovn_lum = ""
      end type input_lum
      type (input_lum) :: in_lum

!! calibration
      type input_chg
       character(len=25) :: cal_parms = ""
       character(len=25) :: cal_upd = ""
       character(len=25) :: codes_sft = ""
       character(len=25) :: wb_parms_sft = ""
       character(len=25) :: water_balance_sft = ""
       character(len=25) :: ch_sed_budget_sft = ""
       character(len=25) :: ch_sed_parms_sft = ""
       character(len=25) :: plant_parms_sft = ""
       character(len=25) :: plant_gro_sft = ""
      end type input_chg
      type (input_chg) :: in_chg

!! initial conditions
      type input_init
       character(len=25) :: plant = ""
       character(len=25) :: soil_plant_ini = ""
       character(len=25) :: om_water = ""
       character(len=25) :: pest_soil = ""
       character(len=25) :: pest_water = ""
       character(len=25) :: path_soil = ""
       character(len=25) :: path_water = ""
       character(len=25) :: hmet_soil = ""
       character(len=25) :: hmet_water = ""
       character(len=25) :: salt_soil = ""
       character(len=25) :: salt_water = ""
      end type input_init
      type (input_init) :: in_init

!! soils
      type input_soils
       character(len=25) :: soils_sol = ""
       character(len=25) :: nut_sol = ""
       character(len=25) :: lte_sol = ""
       character(len=25) :: lyr_depths = ""
      end type input_soils
      type (input_soils) :: in_sol

!! decision_table
      type input_condition
       character(len=25) :: dtbl_lum = ""
       character(len=25) :: dtbl_res = ""
       character(len=25) :: dtbl_scen = ""
       character(len=25) :: dtbl_flo = ""
      end type input_condition
      type (input_condition) :: in_cond

!! regions
      type input_regions
        character(len=25) :: ele_lsu = ""
        character(len=25) :: def_lsu = ""
        character(len=25) :: ele_reg = ""
        character(len=25) :: def_reg = ""
        character(len=25) :: cal_lcu = ""
        character(len=25) :: ele_cha = ""
        character(len=25) :: def_cha = ""
        character(len=25) :: def_cha_reg = ""
        character(len=25) :: ele_aqu = ""
        character(len=25) :: def_aqu = ""
        character(len=25) :: def_aqu_reg = ""
        character(len=25) :: ele_res = ""
        character(len=25) :: def_res = ""
        character(len=25) :: def_res_reg = ""
        character(len=25) :: ele_psc = ""
        character(len=25) :: def_psc = ""
        character(len=25) :: def_psc_reg = ""
      end type input_regions
      type (input_regions) :: in_regs

!! carbon
      type input_carbon
        character(len=25) :: basins_cbn = ""
        character(len=25) :: coefficients_cbn = ""
        character(len=25) :: co2_yr_cbn = ""
      end type input_carbon
      type (input_carbon) :: in_carbon

!! salt
      type input_salt
        character(len=25) :: init_slt = ""
        character(len=25) :: channel_slt = ""
        character(len=25) :: hru_slt = ""
        character(len=25) :: fert_slt = ""
        character(len=25) :: irr_slt = ""
        character(len=25) :: plants_slt = ""
        character(len=25) :: road_slt = ""
        character(len=25) :: urban_slt = ""
        character(len=25) :: uptake_slt = ""
        character(len=25) :: res_slt = ""
      end type input_salt
      type (input_salt) :: in_salt

!! constituents
      type input_constituents
        character(len=25) :: cs_db = ""
        character(len=25) :: init_cs = ""
        character(len=25) :: channel_cs = ""
        character(len=25) :: hru_cs = ""
        character(len=25) :: fert_cs = ""
        character(len=25) :: irr_cs = ""
        character(len=25) :: plants_boron_cs = ""
        character(len=25) :: reactions_cs = ""
        character(len=25) :: uptake_cs = ""
        character(len=25) :: urban_cs = ""
        character(len=25) :: res_cs = ""
        character(len=25) :: streamobs_cs = ""
        character(len=25) :: init_aqu_cs = ""
        character(len=25) :: init_cha_cs = ""
        character(len=25) :: reservoir_cs = ""
        character(len=25) :: wetland_cs = ""
        character(len=25) :: nut_rte = ""
      end type input_constituents
      type (input_constituents) :: in_constit

!! manure
      type input_manure
        character(len=25) :: manure_frt = ""
        character(len=25) :: manure_allo = ""
      end type input_manure
      type (input_manure) :: in_manure

!! water_allocation
      type input_water_allocation
        character(len=25) :: wro = ""
        character(len=25) :: water_pipe = ""
        character(len=25) :: water_tower = ""
        character(len=25) :: water_use = ""
        character(len=25) :: water_treat = ""
        character(len=25) :: om_treat = ""
        character(len=25) :: om_use = ""
        character(len=25) :: outside_src = ""
        character(len=25) :: om_osrc = ""
      end type input_water_allocation
      type (input_water_allocation) :: in_wal

!! update
      type input_update
        character(len=25) :: scen_dtl = ""
      end type input_update
      type (input_update) :: in_upd

!! io_path
      type input_path_pcp
        character(len=80) :: pcp = " "
      end type input_path_pcp
      type (input_path_pcp) :: in_path_pcp

     type input_path_tmp
        character(len=80) :: tmp = " "
      end type input_path_tmp
      type (input_path_tmp) :: in_path_tmp

     type input_path_slr
        character(len=80) :: slr = " "
      end type input_path_slr
      type (input_path_slr) :: in_path_slr

     type input_path_hmd
        character(len=80) :: hmd = " "
      end type input_path_hmd
      type (input_path_hmd) :: in_path_hmd

     type input_path_wnd
        character(len=80) :: wnd = " "
      end type input_path_wnd
      type (input_path_wnd) :: in_path_wnd

    type input_path_pet
        character(len=80) :: peti = " "
      end type input_path_pet
      type (input_path_pet) :: in_path_pet

!! gwflow sub-config (populated from gwflow.aqu; defaults = original hardcoded names)
      type input_gwflow
        !! basic
        character(len=25) :: gw_input = "gwflow.input"
        !! cells
        character(len=25) :: hrucell = "gwflow.hrucell"
        character(len=25) :: lsucell = "gwflow.lsucell"
        character(len=25) :: cellhru = "gwflow.cellhru"
        character(len=25) :: huc12cell = "gwflow.huc12cell"
        character(len=25) :: con = "gwflow.con"
        character(len=25) :: chancells = "gwflow.chancells"
        character(len=25) :: rescells = "gwflow.rescells"
        !! exchange
        character(len=25) :: wetland = "gwflow.wetland"
        character(len=25) :: floodplain = "gwflow.floodplain"
        character(len=25) :: canals = "gwflow.canals"
        !! pumping
        character(len=25) :: pumpex = "gwflow.pumpex"
        character(len=25) :: tiles = "gwflow.tiles"
        !! solutes
        character(len=25) :: solutes = "gwflow.solutes"
        character(len=25) :: solutes_minerals = "gwflow.solutes.minerals"
        character(len=25) :: streamobs = "gwflow.streamobs"
        !! observations
        character(len=25) :: hru_pump_observe = "gwflow.hru_pump_observe"
        character(len=25) :: usgs_head = "usgs_annual_head"
      end type input_gwflow
      type (input_gwflow) :: in_gwf

      contains

      end module input_file_module
