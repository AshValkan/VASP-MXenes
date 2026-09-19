"""
neb.py -- CI-NEB setup for ion diffusion on / inside MXene slabs.

Provides the diffusion barrier E_a that Chapter 3 needs to connect the
electronic structure to the *rate* capability (not just the capacity):

    eps_d up  ->  stronger binding   ->  more capacity
    E_a  down ->  faster transport   ->  better rate performance

Those two pull in opposite directions, and showing that V substitution
improves both is what makes the result interesting.

The climbing-image NEB requires the VTST-enabled VASP binary
(``vasp_vtst``); the module writes plain linear-interpolated images by
default and can switch to the IDPP scheme (improved initial guess) when
``pymatgen`` is installed.
"""

from __future__ import annotations

import os
import shutil
from typing import List, Sequence

import numpy as np

try:
    from pymatgen.core import Structure
    from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver
    _HAVE_PMG = True
except Exception:                                    # pragma: no cover
    _HAVE_PMG = False

__all__ = ["interpolate", "make_neb", "neb_incar", "barrier_from_neb"]


# --------------------------------------------------------------------------- #
#  Image construction
# --------------------------------------------------------------------------- #

def _read_poscar(path: str):
    from dband import read_poscar
    st = read_poscar(path)
    inv = np.linalg.inv(st["lattice"])
    cart = st["positions"] @ st["lattice"] if not st["cartesian"] \
        else st["positions"].copy()
    return st, cart, inv


def _write_poscar(path: str, st, cart: np.ndarray, comment: str,
                  sd: str = "T T T"):
    inv = np.linalg.inv(st["lattice"])
    frac = cart @ inv
    order: List[str] = []
    for e in st["elements"]:
        if e not in order:
            order.append(e)
    lines = [comment, "   1.00000000000000"]
    for v in st["lattice"]:
        lines.append("   " + "  ".join(f"{x:16.10f}" for x in v))
    lines.append("  " + "   ".join(order))
    lines.append("  " + "   ".join(
        str(sum(1 for e in st["elements"] if e == s)) for s in order))
    lines.append("Selective Dynamics")
    # Values below are FRACTIONAL (cart @ inv(lattice)).  Labelling them
    # "Cartesian" while writing Direct coordinates would silently give a
    # wrong structure -- VASP does not warn about this.
    lines.append("Direct")
    for e in order:
        for i, el in enumerate(st["elements"]):
            if el == e:
                lines.append("  " + "  ".join(f"{x:16.12f}" for x in frac[i])
                             + f"   {sd}")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def interpolate(initial: str, final: str, n_images: int = 5,
                method: str = "linear", out: str = "NEB") -> List[str]:
    """Builds the intermediate images.

    ``n_images`` counts the *moving* images; the endpoints 00 and
    ``n_images+1`` are always the supplied structures.

    method
      linear -- straight cartesian interpolation (fine for a hop between two
                neighbouring hollows; can start with atoms too close if the
                path passes near the slab).
      idpp   -- image-dependent pair potential pre-relaxation (pymatgen).
                Strongly recommended for intercalation paths through the MXene
                interlayer, where linear interpolation often drives the ion
                straight through a host atom.
    """
    st_i, cart_i, _ = _read_poscar(initial)
    st_f, cart_f, _ = _read_poscar(final)
    if cart_i.shape != cart_f.shape:
        raise ValueError("initial and final structures have different "
                         f"atom counts: {cart_i.shape} vs {cart_f.shape}")
    if (st_i["elements"] != st_f["elements"]):
        raise ValueError("initial and final structures differ in species order")

    images = [cart_i.copy()]
    if method == "idpp" and _HAVE_PMG:
        s_i = Structure(st_i["lattice"], st_i["elements"], cart_i,
                        coords_are_cartesian=True)
        s_f = Structure(st_f["lattice"], st_f["elements"], cart_f,
                        coords_are_cartesian=True)
        images += [s.cart_coords for s in
                   IDPPSolver.from_endpoints([s_i, s_f], n_images + 1,
                                             ).run().images[1:-1]]
    else:
        if method == "idpp":
            print("  [warn] pymatgen not available -> falling back to linear "
                  "interpolation")
        for k in range(1, n_images + 1):
            images.append(cart_i + (cart_f - cart_i) * k / (n_images + 1))
    images.append(cart_f.copy())

    paths = []
    for k, c in enumerate(images):
        d = os.path.join(out, f"{k:02d}")
        p = os.path.join(d, "POSCAR")
        _write_poscar(p, st_i, c, f"NEB image {k:02d}/{len(images)-1:02d}")
        paths.append(p)
    return paths


# --------------------------------------------------------------------------- #
#  INCAR / job assembly
# --------------------------------------------------------------------------- #

def neb_incar(n_images: int, encut: float = 550, ediff: str = "1E-06",
              ediffg: str = "-0.03", potim: float = 0.0,
              spring: float = -5.0, climb: bool = True,
              max_steps: int = 200, spin: bool = True,
              species: Sequence[str] = (), counts: Sequence[int] = (),
              U: dict | None = None, vdw: int = 11,
              system: str = "NEB") -> str:
    """INCAR for a VTST climbing-image NEB run."""
    from alloy import DEFAULT_U
    U = dict(DEFAULT_U if U is None else U)
    lines = [
        f"SYSTEM = {system}",
        "",
        "# ---------- NEB (VTST) ----------",
        f"IMAGES   = {n_images}",
        f"SPRING   = {spring}",
        f"LCLIMB   = {'.TRUE.' if climb else '.FALSE.'}",
        "ICHAIN   = 0",
        "LNEBCELL = .FALSE.",
        "IOPT     = 3",                 # quick-min (1=QM, 3=FIRE recommended)
        "IBRION   = 3",
        f"POTIM    = {potim}",
        f"NSW      = {max_steps}",
        f"EDIFFG   = {ediffg}",
        "",
        "# ---------- electronic ----------",
        f"ENCUT  = {encut}",
        f"EDIFF  = {ediff}",
        f"ISPIN  = {2 if spin else 1}",
        "ISMEAR = 0",
        "SIGMA  = 0.05",
        "LREAL  = AUTO",
        "ALGO   = Normal",
        "NELMIN = 4",
        "NELM   = 200",
        "PREC   = Accurate",
        "GGA    = PE",
        "",
        "# ---------- Hubbard / vdW ----------",
    ]
    if species and U:
        lines += ["LDAU   = .TRUE.", "LASPH  = .TRUE.", "LDAUTYPE = 2",
                  "LDAUPRINT = 2",
                  "LDAUU  = " + "  ".join(f"{U.get(e, 0.0):.2f}" for e in species),
                  "LDAUJ  = " + "  ".join("0.00" for _ in species),
                  "LDAUL  = " + "  ".join(
                      "2" if e in U else "-1" for e in species)]
    if vdw:
        lines.append(f"IVDW   = {vdw}")
    lines += ["",
              "# ---------- IO ----------",
              "ISTART = 0", "INIWAV = 1", "NWRITE = 2",
              "LWAVE  = .FALSE.", "LCHARG = .FALSE.",
              ""]  # dipole correction dropped: DIPOL ref must be slab center
    return "\n".join(lines) + "\n"


def make_neb(initial: str, final: str, out: str = "NEB", n_images: int = 5,
             method: str = "linear", potcar: str | None = None,
             script: str | None = None, species: Sequence[str] = (),
             counts: Sequence[int] = (), **incar_kw) -> str:
    """Creates a complete NEB working directory.

    VTST reads the endpoint structures from ``00/POSCAR`` and
    ``NN/POSCAR`` and every image runs in its own subdirectory, all sharing
    one INCAR/POTCAR/KPOINTS triplet placed in the parent folder.
    """
    paths = interpolate(initial, final, n_images=n_images, method=method,
                        out=out)

    # Species/counts come from the endpoint POSCAR: the LDAU arrays MUST be
    # written in POSCAR order, and getting that wrong is silent (VASP applies
    # the wrong U without complaining).  Never leave it to the caller.
    from dband import read_poscar
    st = read_poscar(initial)
    species = tuple(st["species"])
    counts = tuple(st["counts"])

    with open(os.path.join(out, "INCAR"), "w") as f:
        f.write(neb_incar(n_images, species=species, counts=counts, **incar_kw))

    if potcar and os.path.exists(potcar):
        shutil.copy(potcar, os.path.join(out, "POTCAR"))
    else:
        print("  [warn] no POTCAR supplied -- copy one into", out)

    # KPOINTS: same mesh as the endpoint relaxations
    src = None
    for cand in (os.path.join(os.path.dirname(initial), "KPOINTS"),
                 os.path.join(os.path.dirname(os.path.dirname(initial)),
                              "opt", "KPOINTS")):
        if os.path.exists(cand):
            src = cand
            break
    if src:
        shutil.copy(src, os.path.join(out, "KPOINTS"))
    else:
        with open(os.path.join(out, "KPOINTS"), "w") as f:
            f.write("K-Points\n0\nGamma\n 4  4  1\n 0  0  0\n")

    if script and os.path.exists(script):
        shutil.copy(script, os.path.join(out, "script"))

    # VTST needs the endpoint POSCARs inside 00/ and NN/ (already written by
    # interpolate) -- nothing else to do.
    return out


# --------------------------------------------------------------------------- #
#  Analysis
# --------------------------------------------------------------------------- #

def barrier_from_neb(neb_dir: str = "NEB",
                     endpoint: str = "00/OUTCAR") -> dict:
    """Extracts E_a and the reaction profile from a finished NEB run.

    Falls back to parsing each image's OUTCAR when VTST's ``neb.dat`` /
    ``nebresults.pl`` output is not present.
    """
    from adsorption import read_energy

    energies = []
    for k in sorted(os.listdir(neb_dir)):
        d = os.path.join(neb_dir, k)
        oc = os.path.join(d, "OUTCAR")
        if k.isdigit() and os.path.exists(oc):
            try:
                energies.append((int(k), read_energy(oc)))
            except Exception:                          # noqa: BLE001
                pass
    if not energies:
        return {"images": 0}

    energies.sort()
    e = np.array([v for _, v in energies])
    e_rel = e - e[0]
    imax = int(np.argmax(e_rel))
    return {"images": len(e),
            "energies_rel_EF": [round(float(x), 4) for x in e_rel],
            "E_a_forward": round(float(e_rel.max()), 4),
            "E_a_reverse": round(float(e_rel.max() - e_rel[-1]), 4),
            "dE_reaction": round(float(e_rel[-1]), 4),
            "TS_image": imax}


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="CI-NEB setup / analysis.")
    ap.add_argument("initial", nargs="?", help="POSCAR of the initial state")
    ap.add_argument("final", nargs="?", help="POSCAR of the final state")
    ap.add_argument("-o", "--out", default="NEB")
    ap.add_argument("-n", "--images", type=int, default=5)
    ap.add_argument("--method", default="linear", choices=["linear", "idpp"])
    ap.add_argument("--potcar", default=None)
    ap.add_argument("--script", default=None)
    ap.add_argument("--encut", type=float, default=550)
    ap.add_argument("--analyze", action="store_true",
                    help="instead: analyse an existing NEB directory")
    a = ap.parse_args()

    if a.analyze:
        import json
        print(json.dumps(barrier_from_neb(a.out), indent=2))
    elif a.initial and a.final:
        d = make_neb(a.initial, a.final, out=a.out, n_images=a.images,
                     method=a.method, potcar=a.potcar, script=a.script,
                     encut=a.encut)
        print(f"NEB directory ready: {d}  "
              f"({a.images} moving images + 2 endpoints)")
        print("  Run with the VTST-enabled binary (vasp_vtst) from", d)
    else:
        ap.error("give two POSCARs, or --analyze with -o NEBNDIR")
