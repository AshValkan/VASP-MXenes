#!/bin/bash
# Adaptive core-block aware queue runner for alkali adsorption jobs.
# Blocks (logical CPUs, one block = 5 physical cores = 10 SMT threads):
#   0-9  10-19  20-29  30-39
# A job is launched only when a whole block is free (no vasp_std process
# bound to any CPU inside the block). This prevents PRRTE's default
# behaviour of cramming every job onto the first 5 cores.
BASE=$HOME/calc/VASP-MXenes-alloy
JOBS_FILE=$BASE/ads_alkali_queue.txt
LOG=$BASE/ads_alkali_queue.log
BLOCKS="0-9 10-19 20-29 30-39"

# block "a-b" busy if any vasp_std proc's affinity mask overlaps it
block_busy() {
    local lo=${1%-*} hi=${1#*-} p m mlo mhi
    for p in $(pgrep -x vasp_std); do
        m=$(awk '/Cpus_allowed_list/{print $2}' /proc/$p/status 2>/dev/null)
        [ -z "$m" ] && continue
        case "$m" in
            *-*) mlo=${m%-*}; mhi=${m#*-} ;;
            *)   mlo=$m;     mhi=$m ;;
        esac
        if [ "$mlo" -le "$hi" ] && [ "$mhi" -ge "$lo" ]; then return 0; fi
    done
    return 1
}

# wait for bulk refs first
while [ ! -f $BASE/.refs_done ]; do sleep 30; done

while read -r job; do
    [ -z "$job" ] && continue
    case_name=$(echo "$job" | awk '{print $1}')
    sub=$(echo "$job" | awk '{print $2}')

    # wait for a free core block
    cpus=""
    while [ -z "$cpus" ]; do
        for b in $BLOCKS; do
            if ! block_busy "$b"; then cpus=$b; break; fi
        done
        [ -z "$cpus" ] && sleep 60
    done

    echo "$(date '+%F %T') launching $case_name/$sub on cpus $cpus" >> $LOG
    bash $BASE/run_direct_mpi.sh "$case_name" 5 "$sub" "$cpus" >> $LOG 2>&1
    sleep 150   # let the job spin up so block scanning sees it
done < $JOBS_FILE
echo "$(date '+%F %T') QUEUE DONE" >> $LOG
