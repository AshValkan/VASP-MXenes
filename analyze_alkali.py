#!/usr/bin/env python3
"""Analyze alkali (Li/Na/K) adsorption: E_ads(y) per cation on Ti3-yVyC2O2.

E_ads(M) = [E(slab+2M) - E(slab) - 2*E_bulk_per_atom(M)] / 2   (eV per M)
"""
import csv, re, os

BASE = 'MXenes_alloy'
CASES = [('Ti3V0C2O2', 0.0), ('Ti2.5V0.5C2O2', 0.5),
         ('Ti2V1C2O2', 1.0), ('Ti1.5V1.5C2O2', 1.5)]
SPECIES = ['Li', 'Na', 'K']

def last_E0(path):
    e0 = None
    with open(path) as f:
        for line in f:
            m = re.search(r'E0=\s*([-\d.E+]+)', line)
            if m:
                e0 = float(m.group(1))
    return e0

def converged(path):
    with open(path) as f:
        return 'reached required accuracy' in f.read()

os.makedirs('results', exist_ok=True)
rows = []
for sp in SPECIES:
    e_bulk = last_E0(f'ref_{sp}/OSZICAR')
    e_atom = e_bulk / 2.0
    print(f'E({sp} bcc) = {e_bulk:.5f} eV, per atom = {e_atom:.5f} eV')
    for c, y in CASES:
        d = f'{BASE}/{c}_s2x2_outer/ABC/HM'
        e_slab = last_E0(f'{d}/opt/OSZICAR')
        e_tot = last_E0(f'{d}/ads_{sp}/OSZICAR')
        if e_tot is None:
            print(f'  {sp} y={y}: ads not started/finished, skip')
            continue
        conv = 'Y' if converged(f'{d}/ads_{sp}/OUTCAR') else 'n'
        e_ads = (e_tot - e_slab - 2 * e_atom) / 2.0
        rows.append((sp, y, c, e_slab, e_tot, e_ads, conv))
        print(f'  {sp} y={y}: E_ads = {e_ads:+.3f} eV  (conv={conv})')

with open('results/adsorption_alkali.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['species', 'y', 'case', 'E_slab_eV', 'E_slab_2M_eV',
                'E_ads_per_M_eV', 'converged'])
    for r in rows:
        w.writerow(r)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 4.2), dpi=200)
    colors = {'Li': '#C0392B', 'Na': '#1E8449', 'K': '#1E4FA8'}
    for sp in SPECIES:
        pts = sorted([(r[1], r[5]) for r in rows if r[0] == sp and r[6] == 'Y'])
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], 'o-',
                    color=colors[sp], lw=2, ms=6, label=sp)
    ax.set_xlabel('V content y in Ti$_{3-y}$V$_y$C$_2$O$_2$')
    ax.set_ylabel('E$_{ads}$(M) / eV')
    ax.set_title('Alkali adsorption energy vs V content (Otop site)')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('results/E_ads_alkali_vs_y.png')
    print('wrote results/adsorption_alkali.csv, results/E_ads_alkali_vs_y.png')
except Exception as e:
    print('plot skipped:', e)
