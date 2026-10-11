"""
adsorption.py -- adsorbate placement and adsorption energetics on MXene slabs.

Supports the quantity the thesis needs in Chapter 3:

    E_ads = E(slab + A) - E(slab) - n * E_ref(A)

for A = H, Li, Na, K on the O-terminated  Ti_(3-y)V_yC_2O_2  surface.

Why this module exists rather than reusing upstream ``structure.py``
-------------------------------------------------------------------
``structure.CONTCAR.addT()`` only appends a *termination* at the hollow sites
of a pristine slab.  Adsorption is different in three ways:

1. The adsorbate may sit **on top of an existing -O** (H binds to the surface
   oxygen and turns it into -OH -- this is the accepted picture for proton
   storage in acidic MXene electrodes), not only in a hollow.
2. On an *alloy* the adsorption sites are no longer equivalent.  Which metal
   sits directly beneath a given surface oxygen changes the site energy, and
   that is precisely what the V substitution is supposed to modify.  This
   module can therefore report -- and let you select -- sites by the element
   underneath.
3. Adsorbing on one face only leaves a net dipole, which shifts the total
   energy and corrupts the work function.  ``both_sides=True`` (the default)
   mirrors the adsorbate onto the opposite face so the slab stays symmetric.

Reference states
----------------
For protons in an aqueous acidic electrolyte the community convention is the
computational-hydrogen-electrode one, i.e.  E_ref(H) = 1/2 E(H2).  For metal
cations the reference is the bulk metal.  Both are wired in; fill
:data:`REFERENCES` with *your own* calculated values (they are
pseudopotential- and functional-dependent, so never hard-code someone
else's number in a thesis).
"""

from __future__ import annotations

import os
from typing import Dict, List, Sequence, Tuple

import numpy as np

try:                                    # reuse upstream if available
    from VASPread import OUTCAR
except Exception:                        # pragma: no cover
    OUTCAR = None

__all__ = ["Slab", "read_energy", "adsorption_energy", "REFERENCES"]

#: Default starting heights above the surface plane, in Angstrom.
#: These are *starting* geometries only -- the optimiser moves them.
DEFAULT_HEIGHT = {
    ("H", "Otop"): 0.98,      # O-H bond length
    ("Li", "hollow"): 1.85,
    ("Na", "hollow"): 2.20,
    ("K", "hollow"): 2.60,
    ("Li", "Otop"): 1.90,
    ("Na", "Otop"): 2.25,
    ("K", "Otop"): 2.65,
    ("Zn", "Otop"): 1.95,     # Zn-O bond length
    ("Mg", "hollow"): 2.00,
    ("Zn", "hollow"): 1.90,
}

#: Reference energies per atom, in eV.  **Placeholders -- compute your own.**
#:   H  : 1/2 * E(H2 molecule)     (computational hydrogen electrode)
#:   Li/Na/K : E per atom of the bulk metal
REFERENCES: Dict[str, float] = {
    "H": None, "Li": None, "Na": None, "K": None, "Mg": None, "Zn": None,
}

#: Reference molecules needed to obtain the above, for documentation.
REFERENCE_SPECIES = {
    "H": "1/2 * E(H2 molecule, in a 10x10x10 A box, gamma only)",
    "Li": "E(Li bcc bulk) / 2 atoms",
    "Na": "E(Na bcc bulk) / 2 atoms",
    "K": "E(K bcc bulk) / 2 atoms",
}


# --------------------------------------------------------------------------- #
#  Slab
# --------------------------------------------------------------------------- #

class Slab:
    """A POSCAR/CONTCAR viewed as a slab, with adsorbate helpers."""

    def __init__(self, path: str, m_elements: Sequence[str] = ("Ti", "V"),
                 o_element: str = "O"):
        from dband import read_poscar

        st = read_poscar(path)
        self.path = path
        self.lattice = st["lattice"]
        self.species = st["species"]
        self.counts = st["counts"]
        self.elements: List[str] = st["elements"]
        self.nAtoms = st["nAtoms"]
        self.m_elements = list(m_elements)
        self.o_element = o_element

        frac = np.array(st["positions"])
        cart = frac @ self.lattice if not st["cartesian"] else frac.copy()
        self.frac, self.cart = frac, cart
        self.z = cart[:, 2]
        self.zmin, self.zmax = self.z.min(), self.z.max()
        self.thickness = self.zmax - self.zmin

    # ------------------------------------------------------------ selection --
    def select(self, element: str | None = None) -> List[int]:
        if element is None:
            return list(range(self.nAtoms))
        return [i for i, e in enumerate(self.elements) if e == element]

    def surface(self, element: str, side: str = "top",
                tol: float = 1.2) -> List[int]:
        """Atoms of ``element`` belonging to one face of the slab.

        ``tol`` (Angstrom) bands atoms into the same surface plane, which is
        needed once the surface has relaxed / buckled.
        """
        idx = self.select(element)
        if not idx:
            return []
        z = self.z[idx]
        ref = z.max() if side == "top" else z.min()
        return [i for i, zz in zip(idx, z) if abs(zz - ref) <= tol]

    def metals(self) -> List[int]:
        out = []
        for e in self.m_elements:
            out += self.select(e)
        return out

    # ------------------------------------------------------- site chemistry --
    def site_metal(self, o_index: int) -> Tuple[str | None, float, float]:
        """Which M atom sits directly beneath (or above) a given surface O.

        Returns ``(element, in_plane_distance, delta_z)``.  Because the oxygen
        of an H_M-terminated MXene sits above the *middle* metal layer, a V
        placed on the inner layer is what actually sits under every adsorption
        site -- worth knowing before you interpret the numbers.
        """
        top = self.z[o_index] > (self.zmin + self.thickness / 2)
        o = self.cart[o_index]
        best = (None, 1e9, 0.0)
        for m in self.metals():
            dz = self.z[m] - self.z[o_index]
            if top and dz > -0.05:
                continue                 # want the metal *below* the top O
            if (not top) and dz < 0.05:
                continue
            d = np.hypot(*(self.cart[m][:2] - o[:2]))
            if d < best[1]:
                best = (self.elements[m], float(d), float(dz))
        return best

    def site_inventory(self) -> List[Dict]:
        """One record per surface oxygen, tagged with the metal underneath.

        This is the table that lets you pick "an H on an O above V" versus
        "an H on an O above Ti" -- the comparison that turns the d-band-centre
        argument into a site-specific one.
        """
        out = []
        for side in ("top", "bottom"):
            for o in self.surface(self.o_element, side):
                el, d, dz = self.site_metal(o)
                out.append({"side": side, "o_index": o, "metal_below": el,
                            "in_plane_d": round(d, 3), "dz": round(dz, 3),
                            "z": round(float(self.z[o]), 3)})
        return out

    def hollow_centres(self, side: str = "top") -> List[np.ndarray]:
        """Centroids of the O-O-O three-fold hollows on one face."""
        os_ = self.surface(self.o_element, side)
        if len(os_) < 3:
            return []
        pts = np.array([self.cart[i][:2] for i in os_])
        zref = self.z[os_].max() if side == "top" else self.z[os_].min()

        centres = []
        a = self.lattice[0][:2]
        b = self.lattice[1][:2]
        n1, n2 = 6, 6
        allp = []
        for p in pts:
            for i in range(-1, 2):
                for j in range(-1, 2):
                    allp.append(p + i * a + j * b)
        allp = np.array(allp)
        for p in pts:
            d = np.linalg.norm(allp - p, axis=1)
            near = allp[np.argsort(d)[:3]]
            c = near.mean(axis=0)
            if np.linalg.norm(c - p) < 0.2:
                continue
            if any(np.linalg.norm(c - k) < 0.15 for k in centres):
                continue
            centres.append(c)
        return [np.array([c[0], c[1], zref]) for c in centres]

    # -------------------------------------------------------------- editing --
    def add_atoms(self, element: str, positions: Sequence[Sequence[float]],
                  sd: str = "T T T"):
        """Appends atoms (cartesian, Angstrom) as a new species block."""
        pos = np.asarray(positions, dtype=float)
        if pos.ndim == 1:
            pos = pos.reshape(1, 3)
        self.cart = np.vstack([self.cart, pos])
        self.elements += [element] * len(pos)
        if element in self.species:
            k = self.species.index(element)
            self.counts[k] += len(pos)
        else:
            self.species.append(element)
            self.counts.append(len(pos))
        self.nAtoms = len(self.elements)
        self.frac = self.cart @ np.linalg.inv(self.lattice)
        self.z = self.cart[:, 2]
        self.zmin, self.zmax = self.z.min(), self.z.max()
        self.thickness = self.zmax - self.zmin
        self._sd = sd
        return self

    def write(self, path: str, comment: str = "adsorbate", sd: str = "T T T"):
        inv = np.linalg.inv(self.lattice)
        order: List[str] = []
        for e in self.elements:
            if e not in order:
                order.append(e)
        lines = [comment, "   1.00000000000000"]
        for v in self.lattice:
            lines.append("   " + "  ".join(f"{x:16.10f}" for x in v))
        lines.append("  " + "   ".join(order))
        lines.append("  " + "   ".join(str(sum(1 for e in self.elements
                                               if e == s)) for s in order))
        lines.append("Selective Dynamics")
        # NOTE: the values written below are *fractional* (cart @ inv(lattice)).
        # Labelling this block "Cartesian" while writing Direct coordinates
        # would silently produce a wrong structure -- VASP does not warn.
        lines.append("Direct")
        for e in order:
            for i, el in enumerate(self.elements):
                if el != e:
                    continue
                f = self.cart[i] @ inv
                lines.append("  " + "  ".join(f"{x:16.12f}" for x in f)
                             + f"   {sd}")
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return path


# --------------------------------------------------------------------------- #
#  Placement
# --------------------------------------------------------------------------- #

def place_adsorbate(slab: Slab, adsorbate: str = "H", site: str = "Otop",
                    o_index: int | None = None, metal_below: str | None = None,
                    side: str = "top", height: float | None = None,
                    both_sides: bool = True) -> Dict:
    """Returns the cartesian positions for one (or a mirrored pair of) adsorbate(s).

    Parameters
    ----------
    site : "Otop" (on top of a surface oxygen -- the H/-OH case) or
           "hollow" (above the O-O-O three-fold site -- the cation case).
    o_index : explicit surface oxygen to use.
    metal_below : instead, pick the first site whose underlying metal is this
        element (``"Ti"`` / ``"V"``).  This is what makes the Ti-site vs V-site
        comparison possible.
    both_sides : mirror onto the opposite face, keeping the slab symmetric and
        the net dipole (and hence the spurious energy shift) at zero.
    """
    if height is None:
        height = DEFAULT_HEIGHT.get((adsorbate, site), 1.5)

    if site == "Otop":
        pool = slab.surface(slab.o_element, side)
        if not pool:
            raise ValueError(f"no {slab.o_element} atoms on the {side} face")
        if o_index is None and metal_below is not None:
            for o in pool:
                if slab.site_metal(o)[0] == metal_below:
                    o_index = o
                    break
            if o_index is None:
                raise ValueError(f"no site with {metal_below} underneath")
        if o_index is None:
            o_index = pool[0]
        base = slab.cart[o_index].copy()
        sign = 1.0 if side == "top" else -1.0
        pos = [base + np.array([0.0, 0.0, sign * height])]

    elif site == "hollow":
        cs = slab.hollow_centres(side)
        if not cs:
            raise ValueError("no three-fold hollow found")
        if metal_below is not None:
            chosen, c = None, None
            for cand in cs:
                d = np.linalg.norm(
                    np.array([slab.cart[i][:2] for i in slab.metals()])
                    - cand[:2], axis=1)
                k = int(np.argmin(d))
                if slab.elements[slab.metals()[k]] == metal_below:
                    chosen, c = k, cand
                    break
            if c is None:
                raise ValueError(f"no hollow above {metal_below}")
        else:
            c = cs[0]
        sign = 1.0 if side == "top" else -1.0
        pos = [c + np.array([0.0, 0.0, sign * height])]

    else:
        raise ValueError(f"site must be 'Otop' or 'hollow', got {site!r}")

    if both_sides:
        # mirror through the mid-plane of the slab
        zmid = 0.5 * (slab.zmin + slab.zmax)
        mirrored = []
        for p in pos:
            q = p.copy()
            q[2] = 2 * zmid - p[2]
            mirrored.append(q)
        pos = pos + mirrored

    return {"positions": [p.tolist() for p in pos], "site": site,
            "o_index": o_index, "side": side,
            "metal_below": metal_below or (slab.site_metal(o_index)[0]
                                           if site == "Otop" and o_index is not None
                                           else None),
            "height": height, "both_sides": both_sides,
            "n_adsorbates": len(pos)}


def make_adsorption_case(slab_path: str, out_path: str, adsorbate: str = "H",
                         **kwargs) -> Dict:
    """Reads a slab, adds the adsorbate(s) and writes the new POSCAR."""
    slab = Slab(slab_path)
    info = place_adsorbate(slab, adsorbate=adsorbate, **kwargs)
    slab.add_atoms(adsorbate, info["positions"])
    slab.write(out_path, comment=f"{adsorbate} on {os.path.basename(slab_path)}")
    info["slab"] = slab_path
    info["out"] = out_path
    return info


# --------------------------------------------------------------------------- #
#  Energetics
# --------------------------------------------------------------------------- #

def read_energy(path: str) -> float:
    """Final E0 (eV) from an OUTCAR, or from a plain number file."""
    if os.path.isdir(path):
        path = os.path.join(path, "OUTCAR")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if OUTCAR is not None:
        try:
            e, _ = OUTCAR(path).getEnergy()
            return float(e)
        except Exception:
            pass
    with open(path, "r") as f:                 # fallback: last float in file
        vals = []
        for line in f:
            t = line.split()
            if len(t) == 1:
                try:
                    vals.append(float(t[0]))
                except ValueError:
                    pass
    if not vals:
        raise ValueError(f"no energy parsed from {path}")
    return vals[-1]


def adsorption_energy(E_slab_ads: float, E_slab: float,
                      E_ref: float | None = None, element: str = "H",
                      n_ads: int = 1) -> float:
    """E_ads per adsorbate (eV).  Negative = favourable.

    ``E_ref`` is the energy of one adsorbate atom in its reference state;
    if omitted it is looked up in :data:`REFERENCES`.
    """
    if E_ref is None:
        E_ref = REFERENCES.get(element)
    if E_ref is None:
        raise ValueError(
            f"No reference energy for {element}.  Compute "
            f"{REFERENCE_SPECIES.get(element, 'the reference state')} and pass "
            f"E_ref=... or set adsorption.REFERENCES['{element}'].")
    return (E_slab_ads - E_slab - n_ads * E_ref) / n_ads


def coverage(slab_n_M: int, n_ads: int, one_side: bool = True) -> float:
    """Adsorbates per surface metal atom (ML), for reporting."""
    denom = slab_n_M / 2 if one_side else slab_n_M
    return n_ads / denom if denom else float("nan")


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse
    import json

    ap = argparse.ArgumentParser(description="Adsorption site tools.")
    ap.add_argument("poscar", help="optimised slab POSCAR/CONTCAR")
    ap.add_argument("--inventory", action="store_true",
                    help="list surface O and the metal beneath each")
    ap.add_argument("--adsorbate", default="H")
    ap.add_argument("--site", default="Otop", choices=["Otop", "hollow"])
    ap.add_argument("--metal-below", default=None,
                    help="choose a site whose underlying metal is this element")
    ap.add_argument("--height", type=float, default=None)
    ap.add_argument("--one-side", action="store_true",
                    help="adsorb only on the top face (dipole!)")
    ap.add_argument("-o", "--out", default="POSCAR_ads")
    a = ap.parse_args()

    slab = Slab(a.poscar)

    if a.inventory or not a.adsorbate:
        inv = slab.site_inventory()
        print(f"{len(inv)} surface O sites:")
        for r in inv:
            print("   " + "  ".join(f"{k}={v}" for k, v in r.items()))
        if not a.adsorbate:
            raise SystemExit(0)

    info = make_adsorption_case(a.poscar, a.out, adsorbate=a.adsorbate,
                               site=a.site, metal_below=a.metal_below,
                               height=a.height, both_sides=not a.one_side)
    print(json.dumps({k: v for k, v in info.items() if k != "slab"}, indent=2))
