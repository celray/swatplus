#!/bin/bash
# Checks the OpenMP version after the source has changed, for example after the
# branch has been rebased onto new upstream commits.
#
# usage: tools/omp/check.sh [-r REF] [-t THREADS] [-w WORK] [-c CMAKE_ARGS] MODEL_DIR [MODEL_DIR ...]
#
#   -r REF         commit whose serial results must be matched (default HEAD~1, the
#                  upstream commit under the OpenMP commit)
#   -t THREADS     thread counts for the OpenMP runs (default "1 4")
#   -w WORK        work directory for builds and runs (default ./omp-check)
#   -c CMAKE_ARGS  extra cmake arguments for all three builds (default -DENABLE_NETCDF=ON)
#                  for ifx: -c "-DENABLE_NETCDF=ON -DCMAKE_Fortran_COMPILER=ifx
#                  -DCMAKE_Fortran_FLAGS=-fp-model=consistent", since the reference
#                  commit does not set that flag
#
# Steps
#   1. omp_generate.py writes the OpenMP source; ompaudit.py --check reports shared
#      writes, file writes and locals in the parallel routines that allowed.txt
#      does not list
#   2. builds REF (serial), the working tree (serial) and the working tree (OpenMP)
#      with the same compiler and flags
#   3. uninit_scan.py lists locals that may keep a value from the previous call
#   4. runs every model with the three builds and compares all output files with
#      the REF run; files that differ only in clock times are not counted
# Exit status 0 when nothing is reported, 1 otherwise.
set -u
ref=HEAD~1; threads="1 4"; work=omp-check; cargs="-DENABLE_NETCDF=ON"
while getopts r:t:w:c: o; do
  case $o in
    r) ref=$OPTARG ;; t) threads=$OPTARG ;; w) work=$OPTARG ;; c) cargs=$OPTARG ;;
    *) sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 2 ;;
  esac
done
shift $((OPTIND - 1))
[ $# -ge 1 ] || { sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

top=$(git rev-parse --show-toplevel) || exit 2
tools=$top/tools/omp
mkdir -p "$work" && work=$(realpath "$work")
export OMP_STACKSIZE=${OMP_STACKSIZE:-64M}
jobs=$(nproc 2>/dev/null || echo 4)
status=0

echo "== 1. static checks"
python3 "$tools/omp_generate.py" "$top/src" "$work/src-omp" || exit 1
python3 "$tools/ompaudit.py" "$work/src-omp" hru_control obj_inflow --check "$tools/allowed.txt" || status=1

echo "== 2. builds (reference $ref = $(git rev-parse --short "$ref"))"
rm -rf "$work/ref-src"; mkdir -p "$work/ref-src"
git -C "$top" archive "$ref" | tar -x -C "$work/ref-src" || exit 1
build() {   # build <name> <source dir> [cmake args]
  local n=$1 s=$2; shift 2
  cmake -S "$s" -B "$work/build-$n" -DCMAKE_BUILD_TYPE=Release -DTAG=omp-check $cargs "$@" \
    > "$work/build-$n.log" 2>&1 && cmake --build "$work/build-$n" -j "$jobs" >> "$work/build-$n.log" 2>&1 \
    || { echo "build $n failed, see $work/build-$n.log"; exit 1; }
  ls -t "$work/build-$n"/swatplus-* | head -1
}
exe_ref=$(build ref "$work/ref-src")
exe_ser=$(build serial "$top")
exe_omp=$(build omp "$top" -DENABLE_OPENMP=ON)
[ -x "$exe_ref" ] && [ -x "$exe_ser" ] && [ -x "$exe_omp" ] || exit 1

echo "== 3. locals that may keep a value between calls"
if grep -qs 'COMPILER_ID "GNU"' "$work"/build-serial/CMakeFiles/*/CMakeFortranCompiler.cmake; then
  python3 "$tools/uninit_scan.py" "$top/src" "$work/build-serial" || status=1
else
  echo "skipped: uninit_scan.py needs a gfortran build (run check.sh with gfortran)"
fi

echo "== 4. runs"
norm() {
  sed -E 's#[0-9]{1,2}:[ 0-9]{1,2}:[ 0-9]{1,2}#hh:mm:ss#g; s#[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}#m/d/yyyy#g' "$1"
}
same() {   # same <dir a> <dir b>: prints the files that differ, returns 1 if any
  local a=$1 b=$2 bad=0 f
  while IFS= read -r f; do
    [ -e "$b/$f" ] || { echo "    only in reference: $f"; bad=1; continue; }
    cmp -s "$a/$f" "$b/$f" && continue
    grep -Iq . "$a/$f" && cmp -s <(norm "$a/$f") <(norm "$b/$f") && continue
    echo "    differs: $f"; bad=1
  done < <(cd "$a" && find . -type f | sort)
  while IFS= read -r f; do
    [ -e "$a/$f" ] || { echo "    not in reference: $f"; bad=1; }
  done < <(cd "$b" && find . -type f | sort)
  return $bad
}
run() {   # run <dir> <exe> [env ...]
  local d=$1 e=$2; shift 2
  rm -rf "$d"; mkdir -p "$d"; cp -r "$model"/. "$d"/
  local s=$(date +%s)
  (cd "$d" && env "$@" "$e" > run.log 2>&1) || echo "    run failed: $d"
  echo $(( $(date +%s) - s ))
}
for model in "$@"; do
  model=$(realpath "$model"); m=$(basename "$model"); r=$work/runs/$m
  echo "-- $m"
  t=$(run "$r/ref" "$exe_ref"); echo "  reference serial  ${t}s"
  t=$(run "$r/serial" "$exe_ser" SWATPLUS_BATCH=on); printf '  serial            %ss  ' "$t"
  if same "$r/ref" "$r/serial" > "$r/serial.diff"; then echo identical; else echo DIFFERENT; cat "$r/serial.diff"; status=1; fi
  for n in $threads; do
    t=$(run "$r/omp$n" "$exe_omp" OMP_NUM_THREADS=$n); printf '  openmp %-3s threads %ss  ' "$n" "$t"
    if same "$r/ref" "$r/omp$n" > "$r/omp$n.diff"; then echo identical; else echo DIFFERENT; cat "$r/omp$n.diff"; status=1; fi
  done
done
[ $status -eq 0 ] && echo "== no problems" || echo "== problems reported above"
exit $status
