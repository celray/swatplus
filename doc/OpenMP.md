# OpenMP parallel runs

SWAT+ can run the HRUs of a day in parallel with OpenMP. The output files are the same as those of the serial build made with the same compiler and flags, file for file and byte for byte, for any number of threads.

## Building and running

```
cmake -B build-omp -DCMAKE_BUILD_TYPE=Release -DENABLE_NETCDF=ON -DENABLE_OPENMP=ON
cmake --build build-omp -j
```

The OpenMP build needs Python 3 (see [Generated source](#generated-source)). The executable name ends in `-omp`. On Linux the gfortran runtimes are linked statically and `libgomp` dynamically, since most distributions ship `libgomp` only as a shared library.

| Environment variable | Effect |
|---|---|
| `OMP_NUM_THREADS` | number of threads (default: all cores) |
| `OMP_STACKSIZE` | stack of each thread; use `64M` or more, since local arrays are on the stack |
| `SWATPLUS_BATCH=off` | run every HRU in the serial command loop, as without OpenMP |

Compare results only between builds from the same compiler. gfortran and ifx results differ in the last digits whether or not OpenMP is used. ifx builds use `-fp-model=consistent` (`/fp:consistent` on Windows): with the default model, compiling with OpenMP changes the rounding of some channel results.

ifx 2026.0 gives every thread the data of the main thread when an allocatable `threadprivate` array is passed to `copyin`. `batch_init` therefore allocates the per-thread arrays in each thread.

## Design

### Batches

The command loop (`command.f90`) runs the objects of a day in the order set by `hyd_connect`. `hru_batch_module` groups consecutive HRU commands into batches, where no HRU of a batch receives flow from another HRU of the same batch. The HRUs of a batch get slots 1, 2, ... in command order. `hru_batch_run` runs `obj_inflow` and `hru_control` for every slot in parallel, then, in command order, the steps the serial loop runs after each HRU: the basin sums (`hru_sums`) and the deposition and outflow hydrograph output. The last slot runs on the main thread, so the module scratch left after a batch is what the serial loop leaves after its last HRU.

`obj_inflow.f90` is the block of `command.f90` that moves the current day of the hydrographs and sums the inflow of an object. It was moved out unchanged so a batch can call it for each slot.

### Ordered gate

Some code in `hru_control` touches state that another HRU of the batch can change the same day. It calls `ordered_gate` first, which waits until every earlier slot has finished, so the code sees the state it sees in the serial loop. The gate is used for:

- water allocation (`wallo_control`) and paddy irrigation (`wet_irrp`), which withdraw from shared sources
- scheduled management operations (`mgt_sched`), which write `mgt_out.txt`, basin crop yields and the shared grazing operation
- decision tables that read or write the state of other objects, draw random numbers or write output files (`dtbl_coupled`)
- saturated buffers, septic HRUs (`sep(isep)%opt` is shared by the HRUs of a septic type), pesticide mixing in tillage, `wet_initial` and the `paddy_test.csv` line

Without OpenMP the slots run in order and the gate returns at once.

### Deferred sums

Totals that every HRU adds to every day would make the HRUs finish in strict order if they were gated. They are kept per HRU and added in command order after the batch: the basin sediment budget (`hru_sums`) and the gwflow cell sums of `gwflow_soil` and `gwflow_wetland` (`gwflow_defer_module`). Each sum receives the same terms in the same order as in the serial loop.

### Decision tables

The tables in `dtbl_lum` and `dtbl_res` are shared by all HRUs that use them, and `conditions` writes its results (`act_hit`) into the table. A table that only reads and writes the current HRU runs on a copy owned by the thread (`dtbl_lum_point`, `dtbl_res_point`). Any other table runs on the shared table after the gate. `dtbl_coupled` makes the choice from the condition variables and action types of each table, and returns `.true.` for anything it does not recognise.

### Per-thread scratch

Module variables that `hru_control` uses as scratch are `!$omp threadprivate` in their modules. Each thread then has a copy of them.

### Leftover values

Some values are read before they are written in a call, so the serial loop uses what the previous HRU (or plant, or call) left. Examples are `enratio` and `etday` (reset in `varinit` only as locals of the same name), `delg_p` in `pl_nut_demand` (set only when the plant has maturity heat units), `c_microb_sed` in `nut_orgnc2` and `gsi_adj` in `et_pot`. Each slot starts with these values marked as not written. Each slot publishes what it wrote when it finishes. Where a value can be read unwritten, the code calls the gate, which fills it in from the closest earlier slot that wrote it. The list is in `leftover_get` and `leftover_put` of `hru_batch_module`. Locals of this kind use a module variable under the local name (`use hru_batch_module, only : delg_p => lo_delg_p`), so the routine body is unchanged.

### Generated source

A local declared with an initialiser (`real :: x = 0.`) has the SAVE attribute, so all threads would share it. For the OpenMP build, CMake runs `tools/omp/omp_generate.py`. It writes a copy of `src` to `<build>/src-omp` in which such locals of the parallel routines are declared without the initialiser and assigned at the start of the routine. Nothing generated is committed, and the serial build compiles `src` as it is. Locals that must keep their value between calls are listed as `saved` in `tools/omp/allowed.txt` and are left alone. They are read only after the gate.

## Maintenance

### After upstream changes

Rebase the OpenMP commit onto the new upstream commit and run the checks:

```
git rebase <upstream commit>
tools/omp/check.sh <model folder> [<model folder> ...]
```

Conflicts can only come from the hand-written changes (about 30 files). The generated part is rebuilt from the new source. `check.sh` runs four steps:

1. `ompaudit.py --check` on the generated source lists every module variable written by the parallel routines that is not indexed by the current HRU, not threadprivate and not in `allowed.txt`, every routine that reads or writes files, and every local that keeps its value between calls.
2. It builds the commit under the OpenMP commit (serial), the working tree (serial) and the working tree (OpenMP).
3. `uninit_scan.py` lists initialised locals that gfortran thinks may be read before they are written in a call.
4. It runs each model with the three builds and compares every output file with the run of the reference build.

Use models that switch on as many output files as possible and that exercise the decision tables, water allocation and gwflow options of interest.

### Reports from check.sh

| Report | What to do |
|---|---|
| `shared` | Per HRU under another index name: add it to `allowed.txt`. Scratch: add it to the module's `threadprivate` list. State other HRUs use the same day: put `call ordered_gate` before the code, or keep it per HRU and add it in `hru_sums`. |
| `io` | Put `call ordered_gate` before the writes, then add the routine to `allowed.txt`. |
| `saved` | Make it automatic (the generator does this when it is not listed), or keep it if it is read only after the gate. |
| `local` | Read the paths. If a call can use the previous call's value, make it a leftover value (see above). Otherwise add it to `allowed.txt` with the reason. |
| run differs | Run with `SWATPLUS_BATCH=off` and with 1 thread to see whether the batch order or the threads cause it, then compare the first differing lines. |

## Files

| File | Contents |
|---|---|
| `src/hru_batch_module.f90` | batches, ordered gate, leftover values, decision table copies, `dtbl_coupled` |
| `src/obj_inflow.f90` | inflow of an object (moved from `command.f90`) |
| `src/gwflow_defer_module.f90` | per-HRU logs of the gwflow cell sums |
| `tools/omp/omp_generate.py` | writes the OpenMP source |
| `tools/omp/ompaudit.py`, `tools/omp/fscan.py` | call tree and hazard scan |
| `tools/omp/uninit_scan.py` | locals that may be read before they are written |
| `tools/omp/allowed.txt` | reviewed hazards, with the reason for each |
| `tools/omp/check.sh` | the checks after source changes |
| `test/omp_compare.sh` | serial against OpenMP runs of one model |
