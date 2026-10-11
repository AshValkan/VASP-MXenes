#!/bin/bash
set -e
cd ~/calc/VASP-MXenes-alloy/MXenes_alloy/Ti3V0C2O2_s2x2_outer/ABC/HM
mkdir -p ads_Li_fB
python3 ~/calc/VASP-MXenes-alloy/make_final_state.py ads_Li ads_Li_fB
cp ads_Li/INCAR ads_Li_fB/INCAR
sed -i -E 's/^IBRION *= *2/IBRION  = 1/; s/^AMIX *= *0\.2/AMIX   = 0.1/; s/^AMIX_MAG *= *0\.2/AMIX_MAG = 0.1/' ads_Li_fB/INCAR
grep -E '^IBRION|^AMIX|^EDIFFG|^NSW' ads_Li_fB/INCAR
cp ads_Li/KPOINTS ads_Li/POTCAR ads_Li_fB/
cd ~/calc/VASP-MXenes-alloy && mkdir -p logs
setsid bash run_direct_mpi.sh Ti3V0C2O2_s2x2_outer 5 ads_Li_fB 0-9 > logs/Ti3V0_ads_Li_fB_launch.log 2>&1 < /dev/null &
sleep 25
echo "procs: $(pgrep -c vasp_std)"
tail -3 Ti3V0C2O2_s2x2_outer/ABC/HM/ads_Li_fB/vasp.out 2>/dev/null | cut -c1-70
