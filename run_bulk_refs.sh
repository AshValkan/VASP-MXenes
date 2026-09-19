#!/bin/bash
# Sequential bulk alkali reference runs (tiny jobs, 5 cores), then raise flag.
BASE=/home/atlas/calc/VASP-MXenes-alloy
export LD_LIBRARY_PATH=/home/atlas/calc/pmix-5.0.1/lib:/usr/local/lib
export OMP_NUM_THREADS=1
for sp in Li Na K; do
    cd $BASE/ref_$sp || continue
    echo "$(date '+%F %T') ref_$sp start" >> $BASE/refs.log
    /usr/local/bin/mpirun -np 5 /home/atlas/calc/vasp_std > vasp.out 2> vasp.err
    echo "$(date '+%F %T') ref_$sp exit=$?" >> $BASE/refs.log
done
touch $BASE/.refs_done
echo "$(date '+%F %T') bulk refs DONE" >> $BASE/refs.log
