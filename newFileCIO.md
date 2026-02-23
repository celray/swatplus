# SWAT+ file.cio Redesign

## Overview

The `file.cio` master I/O control file has been redesigned for clarity, consistency, and completeness. The new format uses section-name-based parsing (order-independent), consistent naming conventions, and brings all formerly-hardcoded filenames under file.cio control. A conversion tool is provided to migrate existing datasets.

This work is on the `cio-redesign` branch.

---

## New Format Summary

The redesigned file.cio has 32 sections (plus a header line). Each line begins with a section keyword followed by space-separated filenames. Section order does not matter  - the Fortran reader dispatches by section name.

**New sections added** (not present in the old format):
- `carbon`  - basins, coefficients, CO2 yearly data
- `salt`  - salt transport initialization, sources, and processes (10 files)
- `constituents`  - constituent databases, initialization, and routing (17 files)
- `manure`  - manure fertilizer and allocation
- `water_allocation`  - water rights, pipes, towers, use, treatment, outside sources (9 files)
- `update`  - scenario detail update file

**Sections renamed:**
- `chg` -> `calibration`
- `water_rights` removed  - its sole active file (`water_allocation.wro` -> `water_allocation.wal`) moved into `water_allocation`

**Path handling:**
- Old format used 5 separate lines (`pcp_path`, `tmp_path`, `slr_path`, `hmd_path`, `wnd_path`)
- New format uses a single `io_path` line with 7 positional tokens: pcp, tmp, slr, hmd, wnd, pet, out (use `null` for unused paths)

**Filename conventions:**
- All dashes replaced with underscores (e.g., `weather-sta.cli` -> `weather_sta.cli`)
- Pattern: `<function>.<module_extension>`

### Complete file.cio Layout

Header: `file.cio: written by SWAT+ rev.62.0.0 - Redesigned with complete file references`

| Line | Section | Token 1 | Token 2 | Token 3 | Token 4 | Token 5 | Token 6 | Token 7 | Token 8 | Token 9 | Token 10 | Token 11 | Token 12 | Token 13 | Token 14 | Token 15 | Token 16 | Token 17 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | simulation | time.sim | print.prt | object.prt | object.cnt | | | | | | | | | | | | | |
| 2 | basin | codes.bsn | parameters.bsn | | | | | | | | | | | | | | | |
| 3 | climate | weather_sta.cli | weather_wgn.cli | pet.cli | pcp.cli | tmp.cli | slr.cli | hmd.cli | wnd.cli | atmodep.cli | atmosalt.cli | atmocs.cli | | | | | | |
| 4 | connect | hru.con | hru_lte.con | rout_unit.con | gwflow.con | aquifer.con | aquifer2d.con | channel.con | reservoir.con | recall.con | exco.con | delratio.con | outlet.con | chandeg.con | | | | |
| 5 | channel | initial.cha | channel.cha | hydrology.cha | sediment.cha | nutrients.cha | channel_lte.cha | hyd_sed_lte.cha | temperature.cha | sed_nut.cha | element.ccu | | | | | | | |
| 6 | reservoir | initial.res | reservoir.res | hydrology.res | sediment.res | nutrients.res | weir.res | wetland.wet | hydrology.wet | res_conds.dat | | | | | | | | |
| 7 | routing_unit | rout_unit.def | rout_unit.ele | rout_unit.rtu | rout_unit.dr | | | | | | | | | | | | | |
| 8 | hru | hru_data.hru | hru_lte.hru | | | | | | | | | | | | | | | |
| 9 | exco | exco.exc | exco_om.exc | exco_pest.exc | exco_path.exc | exco_hmet.exc | exco_salt.exc | | | | | | | | | | | |
| 10 | recall | recall.rec | recall.slt | recall.cs | recall_db.rec | pest.com | | | | | | | | | | | | |
| 11 | dr | delratio.del | dr_om.del | dr_pest.del | dr_path.del | dr_hmet.del | dr_salt.del | | | | | | | | | | | |
| 12 | aquifer | initial.aqu | aquifer.aqu | gwflow.aqu | | | | | | | | | | | | | | |
| 13 | herd | animal.hrd | herd.hrd | ranch.hrd | | | | | | | | | | | | | | |
| 14 | link | chan_surf.lin | aqu_cha.lin | | | | | | | | | | | | | | | |
| 15 | hydrology | hydrology.hyd | topography.hyd | field.fld | | | | | | | | | | | | | | |
| 16 | structural | tiledrain.str | septic.str | filterstrip.str | grassedww.str | bmpuser.str | satbuffer.str | | | | | | | | | | | |
| 17 | hru_parm_db | plants.plt | fertilizer.frt | tillage.til | pesticide.pes | metabolite.pes | pathogens.pth | metals.mtl | salts.slt | urban.urb | septic.sep | snow.sno | | | | | | |
| 18 | ops | harv.ops | graze.ops | irr.ops | chem_app.ops | fire.ops | sweep.ops | puddle.ops | transplant.ops | | | | | | | | | |
| 19 | lum | landuse.lum | management.sch | cntable.lum | cons_practice.lum | ovn_table.lum | | | | | | | | | | | | |
| 20 | calibration | cal_parms.cal | calibration.cal | codes.sft | wb_parms.sft | water_balance.sft | ch_sed_budget.sft | ch_sed_parms.sft | plant_parms.sft | plant_gro.sft | | | | | | | | |
| 21 | init | plant.ini | soil_plant.ini | om_water.ini | pest_hru.ini | pest_water.ini | path_hru.ini | path_water.ini | hmet_hru.ini | hmet_water.ini | salt_hru.ini | salt_water.ini | | | | | | |
| 22 | soils | soils.sol | nutrients.sol | soils_lte.sol | soil_lyr_depths.sol | | | | | | | | | | | | | |
| 23 | decision_table | lum.dtl | res_rel.dtl | scen_lu.dtl | flo_con.dtl | | | | | | | | | | | | | |
| 24 | regions | ls_unit.ele | ls_unit.def | ls_reg.ele | ls_reg.def | ls_cal.reg | ch_catunit.ele | ch_catunit.def | ch_reg.def | aqu_catunit.ele | aqu_catunit.def | aqu_reg.def | res_catunit.ele | res_catunit.def | res_reg.def | rec_catunit.ele | rec_catunit.def | rec_reg.def |
| 25 | carbon | basins.cbn | coefficients.cbn | co2_yr.cbn | | | | | | | | | | | | | | |
| 26 | salt | initial.slt | channel.slt | hru.slt | fertilizer.slt | irrigation.slt | plants.slt | road.slt | urban.slt | uptake.slt | reservoir.slt | | | | | | | |
| 27 | constituents | constituents.cs | initial.cs | channel.cs | hru.cs | fertilizer.cs | irrigation.cs | plants_boron.cs | reactions.cs | uptake.cs | urban.cs | reservoir.cs | streamobs.cs | initial_aqu.cs | initial_cha.cs | reservoir_cs.cs | wetland.cs | nutrients.rte |
| 28 | manure | manure.frt | manure_allo.mnu | | | | | | | | | | | | | | | |
| 29 | water_allocation | water_allocation.wal | water_pipe.wal | water_tower.wal | water_use.wal | water_treat.wal | om_treat.wal | om_use.wal | outside_src.wal | om_osrc.wal | | | | | | | | |
| 30 | update | scen_dtl.upd | | | | | | | | | | | | | | | | |
| 31 | io_path | pcp_path | tmp_path | slr_path | hmd_path | wnd_path | pet_path | out_path | | | | | | | | | | |

Use `null` for any unused token position. Section order in the file does not matter.

---

## Source Code Changes

### Core Changes

**`readcio_read.f90`**  - Complete rewrite of the file.cio reader:
- Changed from positional (fixed line order) parsing to section-name-based dispatch via `select case`
- Section names are case-insensitive (lowercased before matching)
- All sections use a common `parse_line_tokens()` tokenizer  - no more Fortran internal reads
- Token positions match the `cio_feedback_design` column layout exactly
- Detects old-format file.cio (missing `io_path` section) and halts with a clear error message directing the user to run the conversion tool

**`input_file_module.f90`**  - Expanded to cover all file.cio entries:
- All default filenames cleared to empty strings  - nothing is assumed; everything comes from file.cio
- Added new derived types: `input_carbon`, `input_salt`, `input_constituents`, `input_manure`, `input_water_allocation`, `input_update`, `input_gwflow`
- Expanded existing types with new fields to match the redesigned format (e.g., `input_cli` gained `atmo_salt`/`atmo_cs`; `input_rec` gained `recall_slt`/`recall_cs`/`recall_db`/`pest_com`)
- Replaced `input_water_rights`/`in_watrts` with `input_water_allocation`/`in_wal`

**`proc_bsn.f90`**  - Fixed initialization order:
- `readcio_read` now called before any file readers (`basin_read_cc`, `basin_read_objs`, `time_read`)
- Added call to `gwflow_aqu_read` for the gwflow sub-config

**`utils.f90`**  - Added 3-tier `open_file()` utility:
- `open_file(iunit, filename, required, hardcoded)`  - centralised file-open logic
- **Required** files (`required=.true.`): prints error to console + diagnostics.out + simulation.out, then `stop 1`
- **Standard** files (listed in file.cio): prints warning to console + diagnostics.out + simulation.out, returns `.false.`
- **Optional/hardcoded** files: prints note to diagnostics.out + simulation.out only (silent on console), returns `.false.`
- Null or empty filenames are silently skipped (no warning)
- Also added `validate_array_allocation()` for post-read integrity checks

### File Reader Conversions (~116 files)

Every file reader that previously used `inquire(file="hardcoded", exist=i_exist)` + manual `open()` has been converted to use `open_file()`. This covers approximately 116 individual inquire+open patterns across ~112 source files.

The conversion followed one of two patterns:
- **Pattern A** (allocate fallback): `if (open_file(unit, variable)) then ... else allocate(array(0:0)) end if`
- **Pattern B** (no fallback): `if (open_file(unit, variable)) then ... end if`

Key categories of converted files:

| Category | Files | Examples |
|---|---|---|
| Carbon | 3 | `carbon_read.f90`, `carbon_coef_read.f90`, `co2_read.f90` |
| Salt | 8 | `salt_hru_read.f90`, `salt_cha_read.f90`, `salt_irr_read.f90`, etc. |
| Constituents | 9 | `cs_hru_read.f90`, `cs_aqu_read.f90`, `cs_cha_read.f90`, etc. |
| Manure | 2 | `manure_parm_read.f90`, `manure_allocation_read.f90` |
| Water allocation | 6 | `water_use_read.f90`, `water_pipe_read.f90`, etc. |
| Operations | 2 | `mgt_read_puddle.f90`, `plant_transplant_read.f90` |
| Channel & Reservoir | 20 | `ch_read.f90`, `res_read_hyd.f90`, `sd_hydsed_read.f90`, etc. |
| Aquifer & Climate | 12 | `aqu_read.f90`, `cli_staread.f90`, `cli_wgnread.f90`, etc. |
| Database/Parameter | 13 | `plant_parm_read.f90`, `fert_parm_read.f90`, `soil_db_read.f90`, etc. |
| Management/Scenarios | 16 | `mgt_read_harvops.f90`, `scen_read_filtstrip.f90`, `sdr_read.f90`, etc. |
| Calibration/Decision | 13 | `cal_parm_read.f90`, `dtbl_lum_read.f90`, `dtbl_res_read.f90`, etc. |
| Salt/CS/Nutrients | 15 | `salt_fert_read.f90`, `cs_fert_read.f90`, `pest_cha_res_read.f90`, etc. |
| Misc | 14 | `wet_read.f90`, `hru_lte_read.f90`, `om_osrc_read.f90`, etc. |

### Hardcoded Filenames Eliminated

All 18 formerly-hardcoded filenames (files that existed on disk but were never listed in the old file.cio) have been wired to file.cio variables. Six new token slots were added to accommodate files that previously had no file.cio entry:

| Section | New Token | Variable | Source File |
|---|---|---|---|
| recall (pos 4) | `recall_db.rec` | `in_rec%recall_db` | `recall_read.f90` |
| recall (pos 5) | `pest.com` | `in_rec%pest_com` | `recall_read.f90` |
| reservoir (pos 9) | `res_conds.dat` | `in_res%res_conds` | `res_read_conds.f90` |
| soils (pos 4) | `soil_lyr_depths.sol` | `in_sol%lyr_depths` | `soils_init.f90` |
| water_allocation (pos 8) | `outside_src.wal` | `in_wal%outside_src` | `water_osrc_read.f90` |
| water_allocation (pos 9) | `om_osrc.wal` | `in_wal%om_osrc` | `om_osrc_read.f90` |

Zero hardcoded filenames remain in the codebase. Every input file is now controlled through file.cio.

### Climate Station Files  - io_path Support

The 6 climate measurement readers (`cli_pmeas.f90`, `cli_tmeas.f90`, `cli_smeas.f90`, `cli_hmeas.f90`, `cli_wmeas.f90`, `cli_petmeas.f90`) now correctly prepend the io_path prefix to both the station list file and individual timeseries files. Missing station files produce warnings via `open_file()` instead of silently allocating empty arrays.

### gwflow.aqu Sub-Config Reader

A new reader (`gwflow_aqu_read.f90`) parses a `gwflow.aqu` sub-config file referenced from file.cio. This replaces ~30 hardcoded gwflow filenames across 4 source files (`gwflow_read.f90`, `gwflow_chan_read.f90`, `basin_read_objs.f90`, `wet_read_hyd.f90`). The `input_gwflow` type (`in_gwf`) holds 18 configurable filenames covering basic config, cell mapping, exchange, pumping, solutes, and observations. If `gwflow.aqu` is absent, defaults matching the original hardcoded names are used for backward compatibility.

### Required File Enforcement

8 critical file readers now halt cleanly when a required file is listed in file.cio but missing from disk:
- `time_read.f90`, `hru_read.f90`, `recall_read.f90`, `dr_db_read.f90`, `res_read.f90`, `ru_read.f90`, `ru_read_elements.f90`, `aqu2d_read.f90`

Previously these would either crash with a floating-point exception or continue in an undefined state.

### Old-Format Detection

`readcio_read.f90` tracks whether `io_path` was encountered during parsing. If it was not (indicating an old-format file.cio), the model halts with a clear message directing the user to run the conversion tool.

### Other Fixes

- `swift_output.f90`  - Updated `in_sim%cs_db` reference to `in_constit%cs_db` (canonical location)
- `soils_init.f90`  - Hoisted `soil_lyr_depths.sol` open outside the soil-type loop (was producing 260 duplicate diagnostic messages)
- `water_allocation_read.f90`  - Updated to use `in_wal%wro` instead of the removed `in_watrts%transfer_wro`

---

## Conversion Tool

A conversion tool is provided in both Python and C++ to migrate old-format datasets. The tool is maintained in a separate repository: https://github.com/celray/fileCioConverter

**Usage:**
```
convertCIO [path]
```
- `path` can be a TxtInOut directory or a direct file path; defaults to current directory
- Creates a backup as `file.cio.old.format`
- Converts in place and renames dash-filenames on disk
- Scans the data directory for formerly-hardcoded files and adds them to the appropriate sections
- Detects already-converted files and exits cleanly (idempotent)

**Transformations:**
1. `chg` section renamed to `calibration`
2. `water_rights` section removed; `water_allocation.wro` -> `water_allocation.wal` moved to `water_allocation`
3. Old path lines consolidated into single `io_path` line
4. Dashes replaced with underscores in filenames (both in file.cio and on disk)
5. New sections added with appropriate entries (scanned from disk)
6. Token counts adjusted per section for new fields
7. `cs_db` from old simulation line moved to `constituents` section

---

## Testing

| Test | Result |
|---|---|
| 42-year baseline simulation (new format) | Pass |
| Required file missing -> clean halt | Pass |
| Listed file missing -> warning + continue | Pass |
| Null entries -> silently skipped | Pass |
| Conversion tool round-trip (old -> convert -> simulate) | Pass |
| Case-insensitive section names | Pass |
| Section order independence | Pass |
| Conversion tool idempotency | Pass |
| Old-format file.cio detection -> clear error | Pass |

Tests not yet run: io_path with real directories, new sections with real salt/carbon/gwflow data, gwflow.aqu sub-config with a gwflow dataset.

---

## Open Items

- **`water_allocation.wal` naming**: Using `.wal` extension for now; final name TBD pending review with Jeff/Nancy. The filename stem matching the section name may be confusing.
- **Conversion tool token updates**: The Python and C++ tools need updating for the 6 new file.cio token positions (recall 4-5, reservoir 9, soils 4, water_allocation 8-9).
- **Simulation results comparison**: Both swatplus-61 (unmodified) and swatplus-62 (cio-redesign) ran 42-year simulations successfully. Output comparison pending.

---

## Files Changed (src/)

162 source files modified. Key files:

| File | Nature of Change |
|---|---|
| `readcio_read.f90` | Complete rewrite  - section-name dispatch, token parser |
| `input_file_module.f90` | Expanded with new types, cleared all defaults |
| `proc_bsn.f90` | Initialization order fix, gwflow_aqu_read call |
| `utils.f90` | New open_file() and validate_array_allocation() |
| `gwflow_aqu_read.f90` | New file  - gwflow sub-config reader |
| ~157 other `*_read.f90` files | Converted to open_file(), hardcoded filenames wired to file.cio |
