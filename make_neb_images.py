#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build CI-NEB image set (linear interpolation) for Li migration.
Usage: python3 make_neb_images.py <HM_dir>
Creates <HM_dir>/neb_Li/{00..04}/POSCAR + root INCAR/KPOINTS/POTCAR."""
import sys, numpy as np
hm = sys.argv[1]

def load(path):
    L = open(path).read().splitlines()
    s = float(L[1]); A = np.array([[float(x) for x in L[2+i].split()[:3]] for i in range(3)])*s
    sp = L[5].split(); ct = [int(x) for x in L[6].split()]
    mi = next(i for i in range(7,12) if L[i].strip() and L[i].strip()[0].upper() in 'DCS')
    if L[mi].strip()[0].upper()=='S': mi += 1
    assert L[mi].strip()[0].upper()=='D', 'expect Direct coords'
    n = sum(ct)
    D = np.array([[float(x) for x in L[mi+1+i].split()[:3]] for i in range(n)])
    return L[0], s, A, sp, ct, D

h0, s, A, sp, ct, D0 = load(f'{hm}/ads_Li/CONTCAR')
h1, _, _, _, _, D1 = load(f'{hm}/ads_Li_fB/CONTCAR')
# align fractional coords: nearest image per atom (avoid jump across cell)
D1 = D0 + (D1 - D0 - np.round(D1 - D0))
import os
root = f'{hm}/neb_Li'
os.makedirs(root, exist_ok=True)
NIMG = 3
for k in range(NIMG + 2):
    d = f'{root}/{k:02d}'
    os.makedirs(d, exist_ok=True)
    f = k / (NIMG + 1)
    Dk = D0*(1-f) + D1*f
    with open(f'{d}/POSCAR', 'w') as fo:
        fo.write(f'{h0} image {k:02d} (f={f:.3f})\n{s}\n')
        for row in A/s:
            fo.write(' '.join(f'{x:16.10f}' for x in row) + '\n')
        fo.write(' '.join(sp) + '\n' + ' '.join(str(c) for c in ct) + '\nDirect\n')
        for row in Dk:
            fo.write(' '.join(f'{x:16.10f}' for x in row) + '\n')
print(f'images 00..{NIMG+1:02d} written (linear interpolation, frac-diff wrapped)')

# root INCAR: base = ads_Li INCAR with NEB overrides
src = open(f'{hm}/ads_Li/INCAR').read().splitlines()
drop = ('IBRION','EDIFF','EDIFFG','NSW','POTIM','IMAGES','SPRING','LCLIMB','ICHAIN','IOPT','ISYM',
        'LCHARG','LWAVE')
keep = [ln for ln in src if not ln.strip().upper().startswith(drop)]
neb = keep + [
 'IBRION  = 1', 'NFREE   = 2', 'POTIM   = 0.1',
 'EDIFF   = 1e-5', 'EDIFFG  = -0.03', 'NSW     = 200',
 'IMAGES  = 3', 'SPRING  = -5', 'LCLIMB  = .TRUE.',
 'ISYM    = 0', 'LCHARG  = .FALSE.', 'LWAVE   = .FALSE.']
open(f'{root}/INCAR','w').write('\n'.join(neb) + '\n')
import shutil
shutil.copy(f'{hm}/ads_Li/KPOINTS', f'{root}/KPOINTS')
with open(f'{root}/POTCAR','wb') as fo:
    for p in (f'{hm}/ads_Li/POTCAR',):
        fo.write(open(p,'rb').read())
print('root INCAR/KPOINTS/POTCAR written')
