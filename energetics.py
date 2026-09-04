"""
energetics.py -- formation / mixing energies and configurational ranking
for alloyed MXenes.

Everything here works on *total energies from finished VASP runs*, so it can
be used as soon as the optimisations come back -- no DOS required.

The three quantities the thesis needs:

1. **Formation energy** E_form, to show that the alloy is a plausible phase
   at all (and to compare against the competing decomposition into binaries).

2. **Mixing energy** E_mix, the convex-hull ordinate:

       E_mix(y) = E[Ti_(3-y)V_yC_2O_2]
                  - (1 - y/3) * E[Ti_3C_2O_2]
                  - (y/3)     * E[V_3C_2O_2]

   E_mix < 0 means the solid solution is favoured over phase separation.
   Note for the write-up: **V_3C_2T_x cannot be synthesised** (it is only a
   computational end member here), so E_mix is a *theoretical* stability
   measure.  Say so explicitly -- a reviewer who works on MXenes will know.

3. **Site-preference energy** -- the difference between the dopant sitting on
   the outer vs the inner M layer, at fixed composition.  This is the number
   that answers "where does V go?", and it is far more convincing when
   reported as an energy than as a picture.

A note on units: all energies are eV.  "per formula unit" divides by the
number of M(n+1)XnT2 units in the supercell, i.e. by ``sx * sy``.
"""

from __future__ import annotations

import os
from typing import Dict, List, Sequence

__all__ = ["formation_energy", "mixing_energy", "site_preference",
           "collect_energies", "rank_configurations", "hull_table"]


def _e(x) -> float:
    """Accepts a number, a path, or a directory containing an OUTCAR."""
    if isinstance(x, (int, float)):
        return float(x)
    from adsorption import read_energy
    return read_energy(x)


def formation_energy(E_compound: float,
                     references: Dict[str, float],
                     stoichiometry: Dict[str, int],
                     n_formula: int = 1) -> float:
    """E_form per formula unit (eV).

    E_form = [ E(compound) - sum_i n_i * mu_i ] / n_formula

    ``references`` maps element -> chemical potential (eV/atom).  For MXenes
    these are conventionally the elemental ground states (hcp Ti, bcc V,
    graphite C, 1/2 E(O2)), which makes E_form negative for stable
    carbides -- a useful sanity check: if you get a *positive* formation
    energy for Ti3C2O2 something is wrong with the reference set.
    """
    E = _e(E_compound)
    mu = sum(references.get(el, 0.0) * n for el, n in stoichiometry.items())
    return (E - mu) / n_formula


def mixing_energy(E_alloy, E_host_end, E_dopant_end,
                  x_dopant: float) -> float:
    """E_mix per formula unit (eV).  Negative -> solid solution favoured.

    ``x_dopant`` is the dopant *fraction of the M sublattice*, i.e.
    ``x = y / (n+1)``.  For Ti2VC2 (n = 2, y = 1) that is x = 1/3.
    """
    E = _e(E_alloy)
    return E - (1 - x_dopant) * _e(E_host_end) - x_dopant * _e(E_dopant_end)


def site_preference(E_outer, E_inner, n_dopant: int = 1) -> float:
    """dE = E(inner) - E(outer), per dopant atom (eV).

    Positive  -> the dopant prefers the **outer** (surface) layer.
    Negative  -> it prefers the **inner** layer.

    Compare this against kT: at 300 K, kT = 0.026 eV, so anything beyond a
    few tenths of an eV is a decisive preference at room temperature.  The
    samples are sintered at 1400 C, where kT = 0.145 eV -- quote that number
    when you argue the observed site occupancy is thermodynamically
    controlled rather than kinetic.
    """
    return (_e(E_inner) - _e(E_outer)) / max(n_dopant, 1)


def collect_energies(root: str = "MXenes_alloy",
                     outcar: str = "opt/OUTCAR") -> List[Dict]:
    """Walks a generated tree and returns [{label, y, pattern, E}, ...]."""
    from dband import parse_label
    from adsorption import read_energy

    rows = []
    for dirpath, _, files in os.walk(root):
        if "OUTCAR" not in files:
            continue
        p = os.path.join(dirpath, "OUTCAR")
        if not dirpath.replace("\\", "/").endswith("opt"):
            continue
        try:
            E = read_energy(p)
        except Exception:                                  # noqa: BLE001
            continue
        parts = os.path.relpath(dirpath, root).split(os.sep)
        row = parse_label(parts[0] if parts else dirpath)
        row["E"] = E
        row["path"] = dirpath
        rows.append(row)
    rows.sort(key=lambda r: (r.get("y", 0.0), r.get("pattern") or ""))
    return rows


def rank_configurations(rows: Sequence[Dict]) -> List[Dict]:
    """Within each composition, ranks the patterns by energy.

    Returns the same rows with ``dE_min`` (eV above the lowest pattern at the
    same y) and ``rank`` added.
    """
    by_y: Dict[float, List[Dict]] = {}
    for r in rows:
        by_y.setdefault(r.get("y", 0.0), []).append(r)
    out = []
    for y, group in by_y.items():
        e0 = min(r["E"] for r in group)
        group = sorted(group, key=lambda r: r["E"])
        for k, r in enumerate(group):
            r = dict(r)
            r["dE_min"] = round(r["E"] - e0, 4)
            r["rank"] = k + 1
            out.append(r)
    out.sort(key=lambda r: (r.get("y", 0.0), r.get("rank", 0)))
    return out


def hull_table(rows: Sequence[Dict], E_host_end: float, E_dopant_end: float,
               n: int = 2) -> List[Dict]:
    """Adds E_mix to every row (needs the two end-member energies)."""
    out = []
    for r in rows:
        y = r.get("y", 0.0)
        x = y / (n + 1)
        r = dict(r)
        r["x_dopant"] = round(x, 4)
        r["E_mix"] = round(mixing_energy(r["E"], E_host_end, E_dopant_end, x), 4)
        out.append(r)
    return out


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse
    import json

    ap = argparse.ArgumentParser(description="Alloy energetics.")
    ap.add_argument("--root", default="MXenes_alloy")
    ap.add_argument("-n", type=int, default=2)
    ap.add_argument("--host-end", type=float, default=None,
                    help="E of the host end member (Ti3C2O2)")
    ap.add_argument("--dopant-end", type=float, default=None,
                    help="E of the dopant end member (V3C2O2, hypothetical)")
    ap.add_argument("-o", "--out", default="energetics.json")
    a = ap.parse_args()

    rows = collect_energies(a.root)
    if not rows:
        raise SystemExit(f"no opt/OUTCAR found under {a.root}")

    rows = rank_configurations(rows)
    if a.host_end is not None and a.dopant_end is not None:
        rows = hull_table(rows, a.host_end, a.dopant_end, n=a.n)

    with open(a.out, "w") as f:
        json.dump(rows, f, indent=2, default=str)

    hdr = f"{'label':<44s}{'pattern':<10s}{'E (eV)':>16s}{'dE (eV)':>10s}"
    if rows and "E_mix" in rows[0]:
        hdr += f"{'E_mix':>10s}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        line = (f"{str(r.get('label')):<44s}{str(r.get('pattern')):<10s}"
                f"{r['E']:>16.6f}{r.get('dE_min', float('nan')):>10.4f}")
        if "E_mix" in r:
            line += f"{r['E_mix']:>10.4f}"
        print(line)
    print(f"\n-> {a.out}")
