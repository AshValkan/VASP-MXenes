#!/bin/bash
# Launch CI-NEB: 3 images x 5 ranks (IMAGES=3 in INCAR, np = 3*5 = 15).
# Usage: run_neb.sh [neb_root] [np] [cpu_block]
#   neb_root : dir containing 00..04/ + INCAR/KPOINTS/POTCAR
#   np       : total ranks, must equal IMAGES * ranks_per_image (default 15)
#   cpu_block: taskset CPU list (default 0-29; PRRTE needs --bind-to none)
set -e
ROOT=${1:-$HOME/calc/VASP-MXenes-alloy/MXenes_alloy/Ti3V0C2O2_s2x2_outer/ABC/HM/neb_Li}
NP=${2:-15}
CPUS=${3:-0-29}
cd "$ROOT"
export LD_LIBRARY_PATH=$HOME/calc/pmix-5.0.1/lib:/usr/local/lib
export OMP_NUM_THREADS=1
setsid taskset -c "$CPUS" /usr/local/bin/mpirun --bind-to none -np "$NP" $HOME/calc/vasp_std > vasp.out 2>&1 < /dev/null &
sleep 20
echo "procs: $(pgrep -c vasp_std)"
ps -o psr= -C vasp_std | sort -n | uniq -c | head -5
