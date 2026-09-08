"""
analyze_alloy.py -- one command from a finished results tree to the thesis
figures of Chapter 3 / Chapter 5.

It glues together :mod:`dband` (electronic structure), :mod:`energetics`
(stability) and the *experimental* numbers you supply, and produces the two
figures that carry the argument:

  1. ``eps_d_vs_y.png``     d-band centre vs V content, with linear fit + R^2.
  2. ``combined_trend.png`` eps_d, E_ads and measured C_g on one normalised
                            axes set -- the "theory predicts, experiment
                            confirms" figure.

It also writes ``dband_table.csv`` / ``energetics.json``, which are the
numbers to paste into the thesis tables.

Usage
-----
    python analyze_alloy.py --root MXenes_alloy \
        --exp measured.json --ads adsorption.json

``measured.json``  : [{"y": 0.0, "C_g": 320.0}, ...]   (F/g, from your GCD data)
``adsorption.json``: [{"y": 0.0, "E_ads": -1.85}, ...] (eV, from adsorption.py)

Both are optional; the script degrades gracefully and tells you what is
missing rather than crashing.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from typing import Dict, List, Sequence

import numpy as np

from dband import PDOS, scan, plot_trend

__all__ = ["load_series", "merge", "write_csv", "run"]


# --------------------------------------------------------------------------- #

def load_series(path: str | None, key: str) -> Dict[float, float]:
    """Reads a JSON/CSV list of ``{y, <key>}`` records into ``{y: value}``."""
    if not path:
        return {}
    if not os.path.exists(path):
        print(f"  [warn] {path} not found -- skipping the {key} series")
        return {}
    if path.lower().endswith(".csv"):
        out = {}
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                try:
                    out[float(r["y"])] = float(r[key])
                except (KeyError, TypeError, ValueError):
                    pass
        return out
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict):                     # {"0.0": 320, ...}
        return {float(k): float(v) for k, v in data.items()}
    return {float(r["y"]): float(r[key]) for r in data
            if r.get("y") is not None and r.get(key) is not None}


def merge(rows: List[Dict], series: Dict[str, Dict[float, float]]
          ) -> List[Dict]:
    """Attaches the external series onto the rows, matched by y."""
    out = []
    for r in rows:
        r = dict(r)
        y = r.get("y")
        for name, table in series.items():
            if y is not None and y in table:
                r[name] = table[y]
        out.append(r)
    return out


def write_csv(rows: List[Dict], path: str):
    if not rows:
        return
    cols: List[str] = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


# --------------------------------------------------------------------------- #

def run(root: str = "MXenes_alloy", n: int = 2,
        m_elements: Sequence[str] = ("Ti", "V"),
        emin: float = -10.0, emax: float = 0.0,
        exp: str | None = None, ads: str | None = None,
        host_end: float | None = None, dopant_end: float | None = None,
        outdir: str = "results", pattern: str | None = None,
        dpi: int = 300) -> Dict:
    """Full analysis pass.  Returns a dict of everything that was produced."""
    os.makedirs(outdir, exist_ok=True)
    produced: Dict[str, str] = {}

    print(f"[1/4] scanning DOSCARs under {root} ...")
    rows = scan(root, m_elements=m_elements, n=n, emin=emin, emax=emax)
    if pattern:
        rows = [r for r in rows if r.get("pattern") == pattern]
    print(f"      {len(rows)} calculations with a DOSCAR")

    print("[2/4] collecting total energies ...")
    from energetics import collect_energies, rank_configurations, hull_table
    try:
        erows = rank_configurations(collect_energies(root))
        if host_end is not None and dopant_end is not None:
            erows = hull_table(erows, host_end, dopant_end, n=n)
        by_label = {str(r.get("label")): r for r in erows}
        for r in rows:
            lab = str(r.get("label"))
            if lab in by_label:
                r["E_total"] = round(by_label[lab]["E"], 6)
                r["dE_min"] = by_label[lab].get("dE_min")
                if "E_mix" in by_label[lab]:
                    r["E_mix"] = by_label[lab]["E_mix"]
        produced["energetics.csv"] = write_csv(erows,
                                                os.path.join(outdir, "energetics.csv")) or ""
        print(f"      {len(erows)} optimised structures")
    except Exception as e:                                     # noqa: BLE001
        print(f"      [skip] energies: {type(e).__name__}: {e}")

    print("[3/4] merging experimental / adsorption series ...")
    series = {"E_ads(H) [eV]": load_series(ads, "E_ads"),
              "C_g [F/g] (measured)": load_series(exp, "C_g")}
    for name, t in series.items():
        print(f"      {name}: {len(t)} points")
    rows = merge(rows, series)

    produced["dband_table.csv"] = write_csv(
        rows, os.path.join(outdir, "dband_table.csv")) or ""
    with open(os.path.join(outdir, "dband_table.json"), "w") as f:
        json.dump(rows, f, indent=2, default=str)

    print("[4/4] plotting ...")
    for key in ("eps_d(Ti)", "eps_d(V)"):
        # Use "enough rows carry the key", not "every row does": the y = 0
        # end member has no V and y = 3 has no Ti, so an all() test would
        # silently suppress the plot for exactly the interesting cases.
        sub = [r for r in rows if r.get(key) is not None]
        if len(sub) < 2:
            continue
        try:
            # Sanitise the WHOLE key -- splitting on "(" would turn both
            # eps_d(Ti) and eps_d(V) into "eps_d" and they would overwrite
            # each other on disk.
            stem = "".join(ch if ch.isalnum() else "_" for ch in key).strip("_")
            p, r2 = plot_trend(
                sub, key=key,
                out=os.path.join(outdir, f"{stem}_vs_y.png"), dpi=dpi)
            produced[os.path.basename(p)] = p
            print(f"      {key}: R^2 = {r2:.4f}  -> {p}")
        except Exception as e:                                 # noqa: BLE001
            print(f"      [skip] {key}: {e}")

    # combined figure: eps_d + whichever extra series are available
    sub = [r for r in rows if r.get("eps_d(Ti)") is not None]
    if len(sub) >= 2:
        extra = {}
        for name in ("E_ads(H) [eV]", "C_g [F/g] (measured)"):
            vals = [r.get(name, np.nan) for r in sub]
            if sum(np.isfinite(v) for v in vals) >= 2:
                extra[name] = vals
        try:
            p, r2 = plot_trend(
                sub, key="eps_d(Ti)", extra=extra or None,
                out=os.path.join(outdir, "combined_trend.png"), dpi=dpi,
                title=r"$\varepsilon_d$ vs experiment, "
                      r"Ti$_{3-y}$V$_y$C$_2$O$_2$")
            produced["combined_trend.png"] = p
            print(f"      combined: R^2(eps_d) = {r2:.4f}  -> {p}")
        except Exception as e:                                 # noqa: BLE001
            print(f"      [skip] combined: {e}")

    print(f"\nDone. Results in ./{outdir}/")
    for k, v in produced.items():
        if v:
            print(f"   {k}")
    return {"rows": rows, "files": produced}


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Analyse a finished alloy-MXene results tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--root", default="MXenes_alloy")
    ap.add_argument("-n", type=int, default=2)
    ap.add_argument("--emin", type=float, default=-10.0)
    ap.add_argument("--emax", type=float, default=0.0)
    ap.add_argument("--exp", default=None,
                    help="JSON/CSV with measured capacitance: y, C_g")
    ap.add_argument("--ads", default=None,
                    help="JSON/CSV with adsorption energies: y, E_ads")
    ap.add_argument("--host-end", type=float, default=None)
    ap.add_argument("--dopant-end", type=float, default=None)
    ap.add_argument("--pattern", default=None,
                    help="restrict to one occupation pattern")
    ap.add_argument("-o", "--outdir", default="results")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    run(root=a.root, n=a.n, emin=a.emin, emax=a.emax, exp=a.exp, ads=a.ads,
        host_end=a.host_end, dopant_end=a.dopant_end, outdir=a.outdir,
        pattern=a.pattern, dpi=a.dpi)


if __name__ == "__main__":
    main()
