#!/usr/bin/env python3
"""Set up alkali (Li/Na/K) adsorption: bulk refs + 4x3 slab+M cases.

Run from repo root with the .venv interpreter.
Reference: E_ads(M) = [E(slab+2M) - E(slab) - 2*E_bulk_per_atom(M)] / 2
"""
import os
import re
import shutil

os.chdir('/home/atlas/calc/VASP-MXenes-alloy')
from adsorption import make_adsorption_case

CASES = ['Ti3V0C2O2', 'Ti2.5V0.5C2O2', 'Ti2V1C2O2', 'Ti1.5V1.5C2O2']
# (species, PP dir, approximate PBE bcc lattice constant / Angstrom)
ALKALI = [('Li', 'Li', 3.49), ('Na', 'Na', 4.21), ('K', 'K_sv', 5.23)]

# ---------------------------------------------------- bulk bcc references
for sp, pp, a in ALKALI:
    ref = f'ref_{sp}'
    os.makedirs(ref, exist_ok=True)
    with open(f'{ref}/POSCAR', 'w') as f:
        f.write(f'bcc {sp} bulk reference\n   1.00000000000000\n')
        for v in ((a, 0, 0), (0, a, 0), (0, 0, a)):
            f.write('   ' + '  '.join(f'{x:16.10f}' for x in v) + '\n')
        f.write(f'  {sp}\n  2\nDirect\n')
        f.write('  0.0  0.0  0.0\n')
        f.write('  0.5  0.5  0.5\n')
    with open(f'{ref}/INCAR', 'w') as f:
        f.write(f'''SYSTEM = bcc {sp} bulk reference
ENCUT  = 550
EDIFF  = 1e-06
EDIFFG = -0.01
ISPIN  = 1
ISMEAR = 0
SIGMA  = 0.05
LREAL  = .FALSE.
ALGO   = Normal
NELM   = 200
PREC   = Accurate
GGA    = PE
IVDW   = 11
IBRION = 2
ISIF   = 3
NSW    = 30
ISTART = 0
INIWAV = 1
LWAVE  = .FALSE.
LCHARG = .FALSE.
''')
    with open(f'{ref}/KPOINTS', 'w') as f:
        f.write('K-Points\n0\nGamma\n13  13  13\n0  0  0\n')
    shutil.copy(f'PP/{pp}/POTCAR', f'{ref}/POTCAR')
    print(ref, 'ready')

# ---------------------------------------------------- slab + M cases
for sp, pp, _ in ALKALI:
    for c in CASES:
        base = f'MXenes_alloy/{c}_s2x2_outer/ABC/HM'
        out = f'{base}/ads_{sp}'
        os.makedirs(out, exist_ok=True)

        info = make_adsorption_case(f'{base}/opt/CONTCAR', f'{out}/POSCAR',
                                    adsorbate=sp, site='Otop', o_index=24,
                                    side='top', both_sides=True)

        incar = open(f'{base}/opt/INCAR').read()
        incar = incar.replace(f'SYSTEM = {c}_s2x2_outer_opt',
                              f'SYSTEM = {c}_s2x2_outer_ads{sp}')
        # closed-shell cation: zero initial moment for the 2 adatoms
        incar = re.sub(r'(MAGMOM\s*=\s*[^\n]+)', r'\1 2*0', incar, count=1)

        # pad +U blocks to the new species count (H-job lesson)
        n = len(open(f'{out}/POSCAR').read().splitlines()[5].split())
        for key, fill in (('LDAUU', '0.00'), ('LDAUJ', '0.00'), ('LDAUL', '-1')):
            def repl(m, key=key, fill=fill, n=n):
                vals = m.group(1).split()
                return key + ' = ' + ' '.join(vals + [fill] * (n - len(vals)))
            incar = re.sub(key + r'\s*=\s*([^\n]+)', repl, incar, count=1)

        open(f'{out}/INCAR', 'w').write(incar)
        shutil.copy(f'{base}/opt/KPOINTS', f'{out}/KPOINTS')

        with open(f'{out}/POTCAR', 'w') as f:
            f.write(open(f'{base}/opt/POTCAR').read())
            f.write(open(f'PP/{pp}/POTCAR').read())
        print(f'{c} ads_{sp} ready (metal_below={info.get("metal_below")})')
print('ALL READY')
