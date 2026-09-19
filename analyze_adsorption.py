#!/usr/bin/env python3
"""Analyze H adsorption results: E_ads(y) trend on Ti3-yVyC2O2 (Otop site).

E_ads = [E(slab+2H) - E(slab) - 2*E_ref(H)] / 2   (per H atom, eV)
E(slab)  = relaxed E0 from opt/OSZICAR
E(slab+2H) = relaxed E0 from ads_H/OSZICAR (takes the last F= line)
E_ref(H) = 0.5 * E(H2) from ref_H2/OSZICAR
"""
import csv, re, sys, os

BASE = 'MXenes_alloy'
CASES = [('Ti3V0C2O2', 0.0), ('Ti2.5V0.5C2O2', 0.5),
         ('Ti2V1C2O2', 1.0), ('Ti1.5V1.5C2O2', 1.5)]

def last_E0(path):
    """Last F= line of OSZICAR -> E0 (column after E0=)."""
    e0 = None
    with open(path) as f:
        for line in f:
            if 'F=' in line or 'E0=' in line:
                m = re.search(r'E0=\s*([-\d.E+]+)', line)
                if m:
                    e0 = float(m.group(1))
    return e0

def ionic_steps(path):
    n = 0
    with open(path) as f:
        for line in f:
            if 'F=' in line:
                n += 1
    return n

def converged(path):
    with open(path) as f:
        return 'reached required accuracy' in f.read()

E_H2 = last_E0('ref_H2/OSZICAR')
E_refH = E_H2 / 2.0
print(f'E(H2) = {E_H2:.5f} eV, E_ref(H) = {E_refH:.5f} eV\n')

rows = []
print(f"{'y':>4} {'E_slab':>10} {'E_slab+2H':>10} {'E_ads(eV)':>10} {'steps':>6} {'conv':>5}")
for c, y in CASES:
    d = f'{BASE}/{c}_s2x2_outer/ABC/HM'
    e_slab = last_E0(f'{d}/opt/OSZICAR')
    e_ads_tot = last_E0(f'{d}/ads_H/OSZICAR')
    steps = ionic_steps(f'{d}/ads_H/OSZICAR')
    conv = 'Y' if os.path.exists(f'{d}/ads_H/OUTCAR') and converged(f'{d}/ads_H/OUTCAR') else 'n'
    e_ads = (e_ads_tot - e_slab - 2 * E_refH) / 2.0
    rows.append((y, c, e_slab, e_ads_tot, e_ads, steps, conv))
    print(f'{y:>4} {e_slab:>10.3f} {e_ads_tot:>10.3f} {e_ads:>10.3f} {steps:>6} {conv:>5}')

os.makedirs('results', exist_ok=True)
with open('results/adsorption_H.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['y', 'case', 'E_slab_eV', 'E_slab_2H_eV', 'E_ads_per_H_eV',
                'ionic_steps', 'converged'])
    for r in rows:
        w.writerow(r)

# plot
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib.rcParams['font.weight'] = 'bold'
    matplotlib.rcParams['axes.labelweight'] = 'bold'
    matplotlib.rcParams['axes.titleweight'] = 'bold'

    ys = [r[0] for r in rows]
    ea = [r[4] for r in rows]
    ea10 = [10.0 * v for v in ea]          # plot in units of 1e-1 eV
    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=200)   # 1:1 aspect
    ax.plot(ys, ea10, 'o-', color='#1E4FA8', lw=2.5, ms=8,
            mfc='#1E4FA8', mec='white', mew=1.5, zorder=3)

    # per-point label offsets tuned so no label crosses the polyline
    offs = {0.0: (0, 11), 0.5: (0, 11), 1.0: (13, -17), 1.5: (16, 5)}
    for x, v in zip(ys, ea):
        dx, dy = offs.get(x, (0, 10))
        ax.annotate(f'{v:.2f}', (x, 10.0 * v), textcoords='offset points',
                    xytext=(dx, dy), ha='center', fontsize=9.5,
                    fontweight='bold')

    ax.set_xlabel('V content y in Ti$_{3-y}$V$_y$C$_2$O$_2$', fontsize=12)
    ax.set_ylabel('E$_{ads}$(H) ($\\times$10$^{-1}$ eV)', fontsize=12)
    ax.set_title('H adsorption energy vs V content (Otop site)',
                 fontsize=12.5, pad=10)
    ax.grid(alpha=0.3)

    # axes: requested x range/ticks; y in 1e-1 eV units, ticks -8..-3
    ax.set_xlim(-0.1, 1.6)
    ax.set_xticks([0, 0.5, 1.0, 1.5])
    ax.set_ylim(-8.8, -2.6)
    ax.set_yticks([-8, -7, -6, -5, -4, -3])

    # bold tick labels, thicker frame and ticks
    for t in ax.get_xticklabels() + ax.get_yticklabels():
        t.set_fontweight('bold')
        t.set_fontsize(11)
    for s in ax.spines.values():
        s.set_linewidth(2.0)
    ax.tick_params(width=2.0, length=6)

    fig.tight_layout()
    fig.savefig('results/E_ads_vs_y.png')
    print('\nwrote results/adsorption_H.csv, results/E_ads_vs_y.png')
except Exception as e:
    print('\nplot skipped:', e)
