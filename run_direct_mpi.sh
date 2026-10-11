#!/bin/bash
# MPI-parallel launcher for VASP cases on <user>@<server>
# Usage: bash run_direct_mpi.sh <CASE_NAME> [NP] [SUB] [CPUS]
#   CASE_NAME : e.g. Ti3V0C2O2_s2x2_outer
#   NP        : MPI ranks per job (default 5; 4 jobs x 5 ranks = 20 cores)
#   SUB       : calculation subdirectory (default opt; e.g. DOS, ads_H, ads_Li)
#   CPUS      : optional logical-CPU block for this job, e.g. "10-19"
#               (one block = 5 physical cores = 10 SMT threads; keeps jobs
#                on disjoint cores -- PRRTE otherwise crams every job onto
#                the first 5 cores of the machine)
# Uses PMIx 5.0.1 to match the /usr/local OpenMPI 5.0.10 + PRRTE stack.
# Skips relaunch if an opt relaxation already converged.
#
# 2026-09-19 FIX: with CPUS, mpirun now gets --bind-to none. PRRTE's default
# bind-to-core policy IGNORES the taskset cpuset inherited from the wrapper
# and rebinds ranks to machine cores 0-4 (verified: 11 stacked jobs, all
# ranks on CPUs 0-9 at ~16% each). --bind-to none lets ranks keep the
# taskset mask and float inside the block (same as verified live rebinds).

CASE=$1
NP=${2:-5}
SUB=${3:-opt}
CPUS=${4:-}
BASE=$HOME/calc/VASP-MXenes-alloy/MXenes_alloy
DIR=$BASE/$CASE/ABC/HM/$SUB
VASP=$HOME/calc/vasp_std

export LD_LIBRARY_PATH=$HOME/calc/pmix-5.0.1/lib:/usr/local/lib:$LD_LIBRARY_PATH
export PATH=/usr/local/bin:$PATH
export OMP_NUM_THREADS=1

cd "$DIR" || { echo "NO_DIR $DIR"; exit 1; }

if [ "$SUB" = "opt" ] && [ -f OUTCAR ] && grep -q "reached required accuracy" OUTCAR; then
    echo "SKIP_CONVERGED $CASE"
    exit 0
fi

if [ -n "$CPUS" ]; then
    setsid taskset -c "$CPUS" /usr/local/bin/mpirun --bind-to none -np $NP "$VASP" \
        > vasp.out 2> vasp.err </dev/null &
else
    setsid /usr/local/bin/mpirun -np $NP "$VASP" \
        > vasp.out 2> vasp.err </dev/null &
fi
echo "LAUNCHED $CASE/$SUB np=$NP cpus=${CPUS:-any} wrapper_pid=$!"
