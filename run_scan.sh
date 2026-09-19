#!/bin/bash
# Launch the 4 outer-pattern relaxations (composition scan) on the server.
# Each job is fully detached (setsid) so it survives the SSH session.
REPO="$HOME/calc/VASP-MXenes-alloy"
VASP="$HOME/calc/vasp_std"
NP=8
CASES="Ti3V0C2O2_s2x2_outer Ti2.5V0.5C2O2_s2x2_outer Ti2V1C2O2_s2x2_outer Ti1.5V1.5C2O2_s2x2_outer"
for c in $CASES; do
  D="$REPO/MXenes_alloy/$c/ABC/HM/opt"
  if [ -f "$D/OUTCAR" ]; then echo "skip $c (OUTCAR exists)"; continue; fi
  cd "$D" || continue
  echo "$(date) launch $c NP=$NP in $D"
  setsid bash -c "export OMP_NUM_THREADS=1; mpirun -np $NP $VASP > vasp.out 2> vasp.err" >/dev/null 2>&1 &
  sleep 2
done
echo "launched 4 relaxations; monitor with: pgrep -fa vasp_std"
