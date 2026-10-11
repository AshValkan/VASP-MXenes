#!/usr/bin/env python3
"""Set up Q3 adsorption calculations: H2 reference + 4 slab+H cases.

Run from the repo root with the .venv interpreter.
"""
import os
import re
import shutil

os.chdir(os.path.expanduser('~/calc/VASP-MXenes-alloy'))
from adsorption import make_adsorption_case

CASES = ['Ti3V0C2O2', 'Ti2.5V0.5C2O2', 'Ti2V1C2O2', 'Ti1.5V1.5C2O2']

# ---------------------------------------------------------------- H2 reference
ref = 'ref_H2'
os.makedirs(ref, exist_ok=True)
a, d = 10.0, 0.74          # box edge (A), H-H bond (A)
with open(f'{ref}/POSCAR', 'w') as f:
    f.write('H2 molecule reference (CHE)\n   1.00000000000000\n')
    for v in ((a, 0, 0), (0, a, 0), (0, 0, a)):
        f.write('   ' + '  '.join(f'{x:16.10f}' for x in v) + '\n')
    f.write('  H\n  2\nDirect\n')
    z1, z2 = 0.5 - d / 2 / a, 0.5 + d / 2 / a
    f.write(f'  {0.5:16.12f} {0.5:16.12f} {z1:16.12f}\n')
    f.write(f'  {0.5:16.12f} {0.5:16.12f} {z2:16.12f}\n')

with open(f'{ref}/INCAR', 'w') as f:
    f.write('''SYSTEM = H2 reference (CHE)
# electronic settings mirror the slab runs for E_ads cancellation
ENCUT  = 550
EDIFF  = 1e-05
EDIFFG = -0.02
ISPIN  = 2
ISMEAR = 0
SIGMA  = 0.05
LREAL  = .FALSE.
ALGO   = Normal
NELM   = 200
PREC   = Accurate
GGA    = PE
IVDW   = 11
# ionic
IBRION = 2
ISIF   = 2
NSW    = 20
# IO
ISTART = 0
INIWAV = 1
LWAVE  = .FALSE.
LCHARG = .FALSE.
MAGMOM = 2*0.5
''')
with open(f'{ref}/KPOINTS', 'w') as f:
    f.write('K-Points\n0\nGamma\n1  1  1\n0  0  0\n')
shutil.copy('PP/H/POTCAR', f'{ref}/POTCAR')
print('ref_H2 ready')

# ------------------------------------------------------------- slab + H cases
for c in CASES:
    base = f'MXenes_alloy/{c}_s2x2_outer/ABC/HM'
    out = f'{base}/ads_H'
    os.makedirs(out, exist_ok=True)

    info = make_adsorption_case(f'{base}/opt/CONTCAR', f'{out}/POSCAR',
                                adsorbate='H', site='Otop', o_index=24,
                                side='top', both_sides=True)
    print(f"{c}: n_H={info['n_adsorbates']} metal_below={info['metal_below']}")

    incar = open(f'{base}/opt/INCAR').read()
    incar = incar.replace(f'SYSTEM = {c}_s2x2_outer_opt',
                          f'SYSTEM = {c}_s2x2_outer_adsH')
    # one MAGMOM line; append the 2 H initial moments
    incar = re.sub(r'(MAGMOM\s*=\s*[^\n]+)', r'\1 2*0.5', incar, count=1)
    open(f'{out}/INCAR', 'w').write(incar)

    shutil.copy(f'{base}/opt/KPOINTS', f'{out}/KPOINTS')

    with open(f'{out}/POTCAR', 'w') as f:
        f.write(open(f'{base}/opt/POTCAR').read())
        f.write(open('PP/H/POTCAR').read())
    print(f'{c}: ads_H ready')
print('ALL READY')
