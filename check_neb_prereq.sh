#!/bin/bash
echo "procs: $(pgrep -c vasp_std 2>/dev/null || echo 0)"
uptime | grep -o 'load average.*'
cd ~/calc/VASP-MXenes-alloy/MXenes_alloy
echo '=== ads_Li CONTCAR sizes ==='
for c in Ti3V0C2O2 Ti2.5V0.5C2O2 Ti2V1C2O2 Ti1.5V1.5C2O2; do
  f=${c}_s2x2_outer/ABC/HM/ads_Li/CONTCAR
  if [ -s $f ]; then echo "$c OK $(stat -c %s $f) bytes"; else echo "$c MISSING"; fi
done
echo '=== ads_Li INCAR key tags (Ti3V0) ==='
grep -E '^IBRION|^EDIFFG|^EDIFF|^AMIX|^NSW|^ISMEAR|^SIGMA|^LDAUU' Ti3V0C2O2_s2x2_outer/ABC/HM/ads_Li/INCAR
echo '=== neb.py head ==='
head -30 ~/calc/VASP-MXenes-alloy/neb.py 2>/dev/null || echo 'no neb.py'
echo '=== vasp NEB strings ==='
strings $HOME/calc/vasp_std 2>/dev/null | grep -m3 -iE 'climbing|SPRING|nudged'
