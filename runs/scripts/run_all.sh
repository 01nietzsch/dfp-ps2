#!/usr/bin/env bash
# Reproduces every QE run used in the submission.
# Usage (from runs/):  QE_BIN=/path/to/qe/bin NP=4 bash scripts/run_all.sh
# Each job runs in its own directory: SCF with the unmodified BN.in, then a
# non-self-consistent 'bands' run that reads that SCF charge density.
set -euo pipefail
cd "$(dirname "$0")/.."                      # runs/
[ -n "${QE_BIN:-}" ] && export PATH="$QE_BIN:$PATH"
NP=${NP:-4}
export OMP_NUM_THREADS=1
MPI="mpirun --allow-run-as-root --oversubscribe -np $NP"

pw.x -h >/dev/null 2>&1 || true
echo "Using $(command -v pw.x)"

run_job () {   # $1 = directory, $2 = bands-type input (or empty), $3 = run bands.x? (yes/no)
  local d=$1 inp=$2 pp=$3
  rm -rf "$d"; mkdir -p "$d"
  cp BN.in B.pz-vbc.UPF N.pz-vbc.UPF "$d"/
  [ -n "$inp" ] && cp "$inp" "$d"/
  ( cd "$d"
    $MPI pw.x -in BN.in > BN.scf.out 2> /dev/null
    grep -q "convergence has been achieved" BN.scf.out || { echo "SCF NOT CONVERGED in $d"; exit 1; }
    if [ -n "$inp" ]; then
      $MPI pw.x -in "$inp" > "${inp%.in}.out" 2> /dev/null
      grep -q "JOB DONE" "${inp%.in}.out" || { echo "bands run failed in $d"; exit 1; }
      # keep all eigenvalues (pw.x does not print them for >100 k-points)
      cp pwscf.save/data-file-schema.xml "${inp%.in}.xml"
    fi
    if [ "$pp" = yes ]; then
      cp ../bands_pp.in .
      $MPI bands.x -in bands_pp.in > bands_pp.out 2> /dev/null
    fi
    rm -rf pwscf.save pwscf.wfc* pwscf.mix* pwscf.xml 2>/dev/null || true
  )
  echo "done: $d"
}

run_job 01_scf_bands  BN.bands.in       yes   # main band structure (d)
run_job 02_dense      BN.bands_dense.in yes   # 3x denser path (gap-location check)
run_job 03_nearK      BN.bands_nearK.in no    # fine path within 10-20% of K
run_job 04_kz         BN.kz.in          no    # k_z dispersion check for (a)

python3 scripts/analyse.py      # numbers for (b)-(d) + all self-checks -> results.txt
python3 scripts/plot_bands.py   # figure -> bands.pdf / bands.png
