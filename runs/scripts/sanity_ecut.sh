#!/usr/bin/env bash
# NOT used for any number in the submission (those all use ecutwfc = 60 Ry as in BN.in).
# Numerical sanity check only: repeat the near-K and k_z runs at ecutwfc = 120 Ry to show that
# the sub-meV wiggles of the valence band near K and the few-meV k_z "dispersion" of the occupied
# bands at 60 Ry come from the finite plane-wave basis, not from physics.
set -euo pipefail
cd "$(dirname "$0")/.."
[ -n "${QE_BIN:-}" ] && export PATH="$QE_BIN:$PATH"
NP=${NP:-4}; export OMP_NUM_THREADS=1
MPI="mpirun --allow-run-as-root --oversubscribe -np $NP"
d=99_sanity_ecut120; rm -rf $d; mkdir -p $d
cp B.pz-vbc.UPF N.pz-vbc.UPF $d/
for f in BN.in BN.bands_nearK.in BN.kz.in; do
  sed 's/ecutwfc = 60.0/ecutwfc = 120.0/' $f > $d/${f%.in}.ecut120.in
done
cd $d
$MPI pw.x -in BN.ecut120.in > BN.scf.ecut120.out 2>/dev/null
grep -q "convergence has been achieved" BN.scf.ecut120.out
cp -r pwscf.save scf.save
for f in BN.bands_nearK BN.kz; do
  rm -rf pwscf.save; cp -r scf.save pwscf.save
  $MPI pw.x -in $f.ecut120.in > $f.ecut120.out 2>/dev/null
  cp pwscf.save/data-file-schema.xml $f.ecut120.xml
done
rm -rf pwscf.save scf.save pwscf.wfc* pwscf.mix* 2>/dev/null || true
echo done
