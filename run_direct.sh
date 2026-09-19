#!/bin/bash
# Run vasp_std directly (singleton MPI, no mpirun) for one case, fully detached.
REPO="$HOME/calc/VASP-MXenes-alloy"
VASP="$HOME/calc/vasp_std"
D="$REPO/MXenes_alloy/$1/ABC/HM/opt"
cd "$D" || exit 1
[ -f OUTCAR ] && { echo "skip $1 (OUTCAR exists)"; exit 0; }
rm -f vasp.out vasp.err
echo "$(date) launch $1 (direct, no mpirun)"
setsid bash -c "export OMP_NUM_THREADS=8; $VASP > vasp.out 2> vasp.err" </dev/null >/dev/null 2>&1 &
disown
echo "launched pid $!"
