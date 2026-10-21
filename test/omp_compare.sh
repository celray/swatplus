#!/bin/bash
# Runs a model with the serial build and with the OpenMP build at several
# thread counts, and checks that every output file is identical.
#
# usage: omp_compare.sh <serial exe> <openmp exe> <model dir> <work dir> [threads ...]
#
#   threads    thread counts to run the OpenMP build with (default: 1 2 4 8)
#
# Each run gets a fresh copy of the model in <work dir>/<run>. Files that differ
# only in clock times or dates (run.log, simulation.out) are not counted.
# Exit status: 0 if all runs match the serial run, 1 if any differs, 2 on bad usage.
#
# Notes
#   - compare builds from the same compiler and flags; gfortran and ifx results
#     differ in the last digits
#   - the OpenMP build keeps large local arrays on the stack; OMP_STACKSIZE is set
#     to 64M unless it is already set
#   - SWATPLUS_BATCH=off makes the OpenMP build run every hru in the serial loop
set -u
if [ $# -lt 4 ]; then
  sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
ser=$(realpath "$1"); omp=$(realpath "$2"); model=$(realpath "$3"); work=$4
shift 4
threads=${*:-1 2 4 8}
export OMP_STACKSIZE=${OMP_STACKSIZE:-64M}

mkdir -p "$work" || exit 2
work=$(realpath "$work")

run() {   # run <name> <exe> [env ...]
  local name=$1 exe=$2; shift 2
  local d=$work/$name
  rm -rf "$d"; mkdir -p "$d"
  cp -r "$model"/. "$d"/
  local s=$(date +%s)
  (cd "$d" && env "$@" "$exe" > run.log 2>&1)
  local rc=$?
  printf '%-10s rc=%d  %ds\n' "$name" "$rc" "$(( $(date +%s) - s ))"
  return $rc
}

norm() {
  sed -E 's#[0-9]{1,2}:[ 0-9]{1,2}:[ 0-9]{1,2}#hh:mm:ss#g; s#[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}#m/d/yyyy#g' "$1"
}

compare() {   # compare <dir a> <dir b>: prints differing files, returns 1 if any
  local a=$1 b=$2 bad=0 f
  while IFS= read -r f; do
    if [ ! -e "$b/$f" ]; then echo "  only in $(basename $a): $f"; bad=1; continue; fi
    cmp -s "$a/$f" "$b/$f" && continue
    if grep -Iq . "$a/$f" && cmp -s <(norm "$a/$f") <(norm "$b/$f"); then continue; fi
    echo "  differs: $f"; bad=1
  done < <(cd "$a" && find . -type f | sort)
  while IFS= read -r f; do
    [ -e "$a/$f" ] || { echo "  only in $(basename $b): $f"; bad=1; }
  done < <(cd "$b" && find . -type f | sort)
  return $bad
}

run serial "$ser" || exit 1
status=0
for n in $threads; do
  run omp$n "$omp" OMP_NUM_THREADS=$n || { status=1; continue; }
  if compare "$work/serial" "$work/omp$n"; then
    echo "omp$n: identical to serial"
  else
    echo "omp$n: DIFFERENT from serial"
    status=1
  fi
done
exit $status
