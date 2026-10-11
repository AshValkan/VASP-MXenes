#!/bin/bash
# Foreground launcher: detaches the adaptive adsorption queue runner and returns.
BASE=$HOME/calc/VASP-MXenes-alloy
if pgrep -f 'bash run_ads_queue.sh' | grep -v $$ > /dev/null 2>&1; then
    echo "queue runner already alive"
else
    setsid nohup bash $BASE/run_ads_queue.sh > /dev/null 2>&1 < /dev/null &
    echo "queue runner launched"
fi
exit 0
