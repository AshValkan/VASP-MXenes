#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the final-state POSCAR for Li migration NEB:
move the top-surface Li to the nearest-neighbour O-top site (PBC-aware).
Usage: python3 make_final_state.py <ads_Li_dir> <out_dir>
Reads  <ads_Li_dir>/CONTCAR (relaxed initial state)
Writes <out_dir>/POSCAR (final-state starting guess)"""
import sys, numpy as np

src, out = sys.argv[1], sys.argv[2]
lines = open(f'{src}/CONTCAR').read().splitlines()
scale = float(lines[1])
A = np.array([[float(x) for x in lines[2+i].split()[:3]] for i in range(3)]) * scale
species = lines[5].split()
counts = [int(x) for x in lines[6].split()]
n = sum(counts)
# robust: locate the coordinate-mode line ('D'irect/'C'art/'S'elective) after counts
mi = None
for i in range(7, 12):
    t = lines[i].strip()
    if t and t[0].upper() in 'DCS':
        mi = i; break
assert mi is not None, 'coord-mode line not found'
mode = lines[mi].strip()[0].upper()
if lines[mi].strip()[0].upper() == 'S':
    mi += 1
    mode = lines[mi].strip()[0].upper()
assert mode == 'D', f'expect Direct coords, got {mode}'
D = np.array([[float(x) for x in lines[mi+1+i].split()[:3]] for i in range(n)])
C = D @ A                                    # cartesian
lo, hi = [], 0
sp_idx = {}
for s, c in zip(species, counts):
    sp_idx[s] = (hi, hi+c); hi += c
# invariant for slab systems: z of cartesian == z of direct*cz only if A[2]=(0,0,cz); assert
assert abs(A[2,0]) < 1e-8 and abs(A[2,1]) < 1e-8, 'unexpected cell orientation'

def xy_minimg(vec):
    # minimum-image via fractional xy
    f = np.linalg.solve(A[:2,:2].T, vec[:2])
    f -= np.round(f)
    return (f @ A[:2,:2].T)[:2]

iLi0, iLi1 = sp_idx['Li']
zLi = C[iLi0:iLi1,2]
top_li = iLi0 + int(np.argmax(zLi))
iO0, iO1 = sp_idx['O']
zO = C[iO0:iO1,2]
ztop = zO.max()
o_top = [iO0+j for j in range(len(zO)) if zO[j] > ztop-0.6]   # top-surface O group
li = C[top_li]
d = [np.linalg.norm(xy_minimg(C[o]-li)) for o in o_top]
o_below = o_top[int(np.argmin(d))]
dd = [(np.linalg.norm(xy_minimg(C[o]-C[o_below])), o) for o in o_top if o != o_below]
dd.sort()
move, o_nb = dd[0]
print(f'top Li z={li[2]:.3f} A; below O={o_below} (d={min(d):.3f}); '
      f'neighbour O={o_nb}, hop distance={move:.3f} A')
assert 2.0 < move < 4.5, 'hop distance out of expected range for O-top lattice'
new = li.copy()
new[:2] = C[o_nb][:2]
C[top_li] = new
D2 = C @ np.linalg.inv(A)
with open(f'{out}/POSCAR', 'w') as f:
    f.write(lines[0] + '\n' + lines[1] + '\n' +
            '\n'.join(' '.join(f'{x:16.10f}' for x in row) for row in A/scale) + '\n')
    f.write(' '.join(species) + '\n' + ' '.join(str(c) for c in counts) + '\nDirect\n')
    for row in D2:
        f.write(' '.join(f'{x:16.10f}' for x in row) + '\n')
print(f'final-state POSCAR -> {out}/POSCAR (Li {top_li}: '
      f'{li[0]:.2f},{li[1]:.2f} -> {new[0]:.2f},{new[1]:.2f}, z kept)')
