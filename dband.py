"""
dband.py -- d-band centre (and higher moments) from VASP DOSCAR files.

This is the module that produces the *central descriptor* of the thesis

        合金化 MXene 电子结构调控与超级电容器储能性能研究

The claim being tested is:

    "Introducing V moves the Ti-site d-band centre (eps_d) toward the Fermi
     level, which strengthens the chemical affinity of the MXene surface for
     electrolyte cations and surface terminations -- and that is the
     electronic-structure origin of the enhanced pseudocapacitance."

For that claim to be more than a picture, eps_d has to be a *number*,
reproducible and comparable across compositions.  This module provides:

  * ``PDOS``          -- element- and site-resolved projected DOS from a
                         DOSCAR + POSCAR pair.  Unlike upstream ``DOS.py`` it
                         does **not** infer the atom list from the SYSTEM tag
                         (which breaks for alloys and supercells); it reads
                         the species/counts lines of the POSCAR, which is
                         authoritative.
  * ``dband_center()``-- eps_d = int(E * rho_d) dE / int(rho_d) dE, with a
                         configurable integration window, per element and per
                         M layer (outer vs inner).
  * ``d_moments()``   -- band width and skewness (2nd/3rd central moments).
                         Skewness matters here: for a fixed band width a
                         positively skewed d band has more weight just below
                         E_F, which is exactly the regime relevant for
                         adsorption.
  * ``scan()``        -- walks a generated results tree and returns the
                         eps_d vs composition table.
  * ``plot_trend()``  -- the thesis figure: eps_d (and optionally E_ads and
                         measured C_g) against V content, with a linear fit
                         and R^2.  Without the R^2 the whole argument stays
                         qualitative and reviewers will say so.

Conventions
-----------
* Energies are shifted so that the Fermi level is at 0 eV.
* ``LORBIT = 11`` (lm-decomposed) is required; the 5 d columns are summed,
  so the exact VASP sub-ordering of the d orbitals is irrelevant here.
* The default window integrates the **occupied** part of the d band only
  (-10 eV .. E_F), which is the standard definition of the d-band centre.
  Set ``emax`` wider if you want the full-band (first-moment) convention --
  but be consistent across every composition, and state which one you used.
"""

from __future__ import annotations

import os
import re
from typing import Dict, Iterable, List, Sequence

import numpy as np

__all__ = ["PDOS", "dband_center", "scan", "plot_trend"]


# --------------------------------------------------------------------------- #
#  POSCAR reader
# --------------------------------------------------------------------------- #

def read_poscar(path: str) -> Dict:
    """Minimal POSCAR/CONTCAR reader -> lattice, species, counts, positions."""
    with open(path, "r") as f:
        raw = f.read().splitlines()

    scale = float(raw[1].split()[0])
    lattice = np.array([[float(x) for x in raw[2 + i].split()[:3]]
                        for i in range(3)]) * scale
    species = raw[5].split()
    counts = [int(x) for x in raw[6].split()]

    i = 7
    mode = raw[i].strip().lower()
    if mode.startswith("s"):                     # Selective dynamics
        i += 1
        mode = raw[i].strip().lower()
    if mode.startswith(("d", "c", "k")):         # Direct / Cartesian
        cart = mode.startswith(("c", "k"))
        i += 1
    else:
        cart = False

    pos, flags = [], []
    for line in raw[i:i + sum(counts)]:
        t = line.split()
        if len(t) < 3:
            break
        pos.append([float(x) for x in t[:3]])
        if len(t) >= 6:
            flags.append(t[3:6])

    pos = np.array(pos, dtype=float)
    elements: List[str] = []
    for s, c in zip(species, counts):
        elements += [s] * c

    return {"lattice": lattice, "species": species, "counts": counts,
            "elements": elements, "positions": pos,
            "cartesian": cart, "flags": flags, "nAtoms": sum(counts)}


# --------------------------------------------------------------------------- #
#  DOSCAR / PDOS
# --------------------------------------------------------------------------- #

class PDOS:
    """Projected density of states, element- and site-resolved.

    Parameters
    ----------
    doscar : path to the DOSCAR file.
    poscar : path to the POSCAR/CONTCAR that generated it.  If omitted it is
        searched for in the DOSCAR folder, then in ``opt/``, then in parent
        directories (see :meth:`_find_poscar`).
    """

    #: VASP lm-decomposed column order (LORBIT = 11), non-spin case.
    ORBITALS_NSP = ("s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "dx2y2")
    D_SLICE = slice(4, 9)          # the five d columns
    P_SLICE = slice(1, 4)

    def __init__(self, doscar: str, poscar: str | None = None):
        self.doscar_path = doscar
        self.folder = os.path.dirname(doscar) or "."

        if poscar is None:
            poscar = self._find_poscar()
        if poscar is None or not os.path.exists(poscar):
            raise FileNotFoundError(
                f"No POSCAR/CONTCAR found in {self.folder} or any parent. "
                f"The atom list cannot be inferred from the SYSTEM tag for "
                f"alloys/supercells -- copy the POSCAR next to the DOSCAR.")

        self.structure = read_poscar(poscar)
        self.elements = self.structure["elements"]
        self.nAtoms = len(self.elements)

        self._read(doscar)

    def _find_poscar(self, max_up: int = 5) -> str | None:
        """Looks for POSCAR/CONTCAR in the DOSCAR folder, then in parents.

        In the generated tree the slab geometry lives in ``opt/POSCAR`` (or as
        ``CONTCAR`` at the case root) while the DOSCAR sits in ``DOS/``, so a
        straight sibling lookup would fail.

        Levels 0-2 of the walk are guaranteed to stay inside one case
        directory, so probing ``opt/`` there can never pick up a different
        composition's geometry; above that only the plain names are checked.
        """
        names = ("POSCAR", "CONTCAR", "poscar", "contcar")
        d = os.path.abspath(self.folder)
        for level in range(max_up + 1):
            for cand in names:
                p = os.path.join(d, cand)
                if os.path.exists(p):
                    return p
            if level <= 2:                       # still inside the same case
                for cand in names:
                    p = os.path.join(d, "opt", cand)
                    if os.path.exists(p):
                        return p
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent
        return None

    # ------------------------------------------------------------------ io --
    def _read(self, path: str):
        with open(path, "r") as f:
            lines = [ln.split() for ln in f]

        # --- locate the "emax emin NEDOS efermi" line ----------------------
        nedos_idx = None
        for i, t in enumerate(lines[:12]):
            # VASP >= 6.x writes 5 fields (emax emin NEDOS efermi weight);
            # older builds write 4.  Accept both.
            if 4 <= len(t) <= 5:
                try:
                    if int(float(t[2])) > 10:
                        nedos_idx = i
                        break
                except ValueError:
                    pass
        if nedos_idx is None:
            raise ValueError(f"Cannot locate the NEDOS line in {path}")

        self.NEDOS = int(float(lines[nedos_idx][2]))
        self.efermi = float(lines[nedos_idx][3])

        # --- where does the total-DOS block start? --------------------------
        # Real VASP DOSCAR: data begins on the line right after the
        # emax/emin/NEDOS/efermi header.  Synthetic/mock files (and a few
        # third-party writers) insert a column-description line such as
        # "energy dos int" -- detect it via a non-numeric first token.
        start = nedos_idx + 1
        if start < len(lines):
            try:
                float(lines[start][0])
            except (ValueError, IndexError):
                start += 1                       # column-description line
        ncol = len(lines[start])
        self.spin = ncol >= 5                    # sp: E, up, dn, iup, idn
        self.nspin_cols = ncol

        stride = self.NEDOS + 1                  # +1 -> separator line

        def block(iat: int) -> np.ndarray:
            a = start + iat * stride
            return np.array([[float(x) for x in r] for r in
                             lines[a:a + self.NEDOS]], dtype=float)

        tot = block(0)
        self.E = tot[:, 0] - self.efermi
        if self.spin:
            self.total = tot[:, 1] + tot[:, 2]
            self.total_up, self.total_dn = tot[:, 1], tot[:, 2]
        else:
            self.total = tot[:, 1]
            self.total_up = self.total_dn = None

        # --- per-atom lm-decomposed blocks ---------------------------------
        # Between blocks VASP >= 6 repeats the emax/emin/NEDOS/efermi header
        # line; some writers emit a blank line instead; some emit nothing.
        # A separator is therefore: an empty line, or a 4-5 token line whose
        # third field round-trips to NEDOS.
        def _is_sep(t) -> bool:
            if not t:
                return True
            if 4 <= len(t) <= 5:
                try:
                    return int(float(t[2])) == self.NEDOS
                except (ValueError, IndexError):
                    return False
            return False

        stride = self.NEDOS + 1
        probe = start + self.NEDOS
        if probe < len(lines) and not _is_sep(lines[probe]):
            stride = self.NEDOS                  # blocks are back-to-back
        self.stride = stride

        norb = 18 if self.spin else 9
        per_atom = np.zeros((self.nAtoms, norb, self.NEDOS), dtype=float)
        for iat in range(self.nAtoms):
            b = block(iat + 1)
            per_atom[iat] = b[:, 1:1 + norb].T
        self.per_atom = per_atom

    # ------------------------------------------------------------ accessors --
    def indices(self, element: str | None = None) -> List[int]:
        """Atom indices belonging to ``element`` (None -> all)."""
        if element is None:
            return list(range(self.nAtoms))
        return [i for i, e in enumerate(self.elements) if e == element]

    def layers(self, n: int = 2, m_elements: Sequence[str] = ("Ti", "V")
               ) -> Dict[str, List[int]]:
        """Splits the M sublattice into its (n+1) layers by cartesian z.

        For n = 2 this returns ``{"outer_bottom", "inner", "outer_top"}``.
        Resolving eps_d per layer is what answers *where* the dopant acts --
        the outer layers are the ones that actually see the electrolyte.
        """
        st = self.structure
        pos = st["positions"]
        if st["cartesian"]:
            z = pos[:, 2]
        else:
            # fractional -> cartesian; column 2 of the lattice is the c axis
            z = (pos @ st["lattice"])[:, 2]
        z = z - z.min()

        m_idx = [i for i, e in enumerate(self.elements) if e in m_elements]
        if not m_idx:
            return {}
        zm = z[m_idx]
        # k-means-free 1-D split: sort and cut at the largest (n) gaps
        order = np.argsort(zm)
        zs = zm[order]
        cuts = np.argsort(np.diff(zs))[-(n):][::-1]
        cuts = np.sort(cuts)
        bounds = [0] + [int(c) + 1 for c in cuts] + [len(zs)]
        names = (["outer_bottom"] + [f"inner{i+1}" for i in range(n - 1)]
                 + ["outer_top"]) if n >= 2 else ["outer_bottom", "outer_top"]

        out: Dict[str, List[int]] = {}
        for k in range(len(bounds) - 1):
            grp = order[bounds[k]:bounds[k + 1]]
            out[names[k] if k < len(names) else f"layer{k}"] = \
                [m_idx[j] for j in grp]
        return out

    def pdos(self, element: str | None = None, orbital: str = "d",
             indices: Sequence[int] | None = None) -> np.ndarray:
        """Summed PDOS over the selected atoms and orbital channel.

        ``orbital`` : "s", "p", "d" or "all".
        """
        if indices is None:
            indices = self.indices(element)
        if not indices:
            return np.zeros(self.NEDOS)

        if self.spin:
            # columns come in (up, dn) pairs per orbital
            blocks = np.zeros((9, self.NEDOS))
            for o in range(9):
                blocks[o] = (self.per_atom[indices, 2 * o, :].sum(axis=0)
                             + self.per_atom[indices, 2 * o + 1, :].sum(axis=0))
        else:
            blocks = self.per_atom[indices, :, :].sum(axis=0)

        if orbital == "s":
            return blocks[0]
        if orbital == "p":
            return blocks[self.P_SLICE].sum(axis=0)
        if orbital == "d":
            return blocks[self.D_SLICE].sum(axis=0)
        if orbital == "all":
            return blocks.sum(axis=0)
        raise ValueError(f"orbital must be s/p/d/all, got {orbital!r}")

    # -------------------------------------------------------- d-band moments --
    def dband_center(self, element: str | None = None,
                     indices: Sequence[int] | None = None,
                     emin: float = -10.0, emax: float = 0.0,
                     orbital: str = "d", norm: bool = True) -> float:
        """eps_d = <E> over the selected d-PDOS inside [emin, emax].

        Parameters
        ----------
        emin, emax : integration window **relative to E_F** (E_F = 0).
            Default (-10, 0) integrates the occupied d band, which is the
            standard d-band-centre definition.
        norm : normalise by the number of selected atoms (so that eps_d can be
            compared between a 2x2 cell with 8 Ti and one with 4 Ti).
        """
        rho = self.pdos(element=element, orbital=orbital, indices=indices)
        E = self.E
        win = (E >= emin) & (E <= emax)
        if not win.any():
            raise ValueError(f"Empty d-band window [{emin}, {emax}]")
        num = np.trapezoid(E[win] * rho[win], E[win])
        den = np.trapezoid(rho[win], E[win])
        if abs(den) < 1e-12:
            return float("nan")
        centre = num / den
        if norm:
            n = len(indices) if indices is not None else len(self.indices(element))
            centre = centre   # centre is already intensive (a ratio); no scaling
        return float(centre)

    def d_moments(self, element: str | None = None,
                  indices: Sequence[int] | None = None,
                  emin: float = -10.0, emax: float = 0.0) -> Dict[str, float]:
        """1st/2nd/3rd moments of the d band: centre, width (RMS), skewness.

        Skewness is worth reporting for this project: at a *fixed* band width,
        a more positively skewed d band places more states immediately below
        E_F, which correlates with stronger adsorbate binding.  If eps_d alone
        fails to correlate with the measured capacitance, skewness is the
        first thing to check.
        """
        rho = self.pdos(element=element, orbital="d", indices=indices)
        E = self.E
        win = (E >= emin) & (E <= emax)
        Ew, rw = E[win], rho[win]
        den = np.trapezoid(rw, Ew)
        if abs(den) < 1e-12:
            return {"center": float("nan"), "width": float("nan"),
                    "skewness": float("nan"), "filling": 0.0}
        c = np.trapezoid(Ew * rw, Ew) / den
        var = np.trapezoid((Ew - c) ** 2 * rw, Ew) / den
        width = float(np.sqrt(max(var, 0.0)))
        m3 = np.trapezoid((Ew - c) ** 3 * rw, Ew) / den
        skew = float(m3 / var ** 1.5) if var > 1e-12 else float("nan")
        return {"center": float(c), "width": width, "skewness": skew,
                "filling": float(den)}

    def dos_at_fermi(self) -> float:
        """N(E_F) from the total DOS -- a proxy for the electronic conductivity."""
        return float(np.interp(0.0, self.E, self.total))

    # ---------------------------------------------------------------- report --
    def summary(self, m_elements: Sequence[str] = ("Ti", "V"), n: int = 2,
                emin: float = -10.0, emax: float = 0.0) -> Dict:
        """Everything the thesis needs from one calculation, in one dict."""
        out = {"folder": os.path.basename(os.path.abspath(self.folder)),
               "nAtoms": self.nAtoms,
               "spin": "sp" if self.spin else "nsp",
               "N(E_F)": round(self.dos_at_fermi(), 4)}
        for el in dict.fromkeys(self.elements):
            if el in m_elements:
                out[f"eps_d({el})"] = round(
                    self.dband_center(el, emin=emin, emax=emax), 4)
                m = self.d_moments(el, emin=emin, emax=emax)
                out[f"w_d({el})"] = round(m["width"], 4)
                out[f"skew_d({el})"] = round(m["skewness"], 4)
        for name, idx in self.layers(n=n, m_elements=m_elements).items():
            if not idx:
                continue
            els = {self.elements[i] for i in idx}
            tag = "/".join(sorted(els)) if len(els) == 1 else "mix"
            out[f"eps_d({tag}@{name})"] = round(
                self.dband_center(indices=idx, emin=emin, emax=emax), 4)
        return out

    # ------------------------------------------------------------------ plot --
    def plot(self, elements: Iterable[str] | None = None,
             out: str | None = None, orbital: str = "d",
             xlim: Sequence[float] = (-10, 6),
             show_center: bool = True, dpi: int = 300,
             emin: float = -10.0, emax: float = 0.0):
        """PDOS plot with the d-band centre marked for each element."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if elements is None:
            elements = [e for e in dict.fromkeys(self.elements)]
        colors = plt.cm.tab10(np.linspace(0, 1, 10))

        fig, ax = plt.subplots(figsize=(5.2, 3.4))
        ax.plot(self.E, self.total, color="0.75", lw=1.0, label="Total")

        for k, el in enumerate(elements):
            rho = self.pdos(element=el, orbital=orbital)
            ax.plot(self.E, rho, color=colors[k % 10], lw=1.3,
                    label=f"{el}-{orbital}")
            if show_center and el in ("Ti", "V", "Nb", "Mo", "Cr"):
                c = self.dband_center(el, emin=emin, emax=emax)
                ax.axvline(c, color=colors[k % 10], ls="--", lw=1.0)
                ax.annotate(f"$\\varepsilon_d$={c:+.3f}", xy=(c, ax.get_ylim()[1]),
                            xytext=(3, -4 - 9 * k), textcoords="offset points",
                            fontsize=7, color=colors[k % 10])

        ax.axvline(0, color="k", lw=0.9, ls=":")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlim(*xlim)
        ax.set_xlabel(r"$E - E_F$ (eV)")
        ax.set_ylabel(f"{orbital}-PDOS (states/eV)")
        ax.legend(frameon=False, fontsize="x-small")
        fig.tight_layout()

        if out is None:
            out = os.path.join(self.folder, f"pdos_{orbital}.png")
        fig.savefig(out, dpi=dpi)
        plt.close(fig)
        return out


# --------------------------------------------------------------------------- #
#  Functional alias
# --------------------------------------------------------------------------- #

def dband_center(doscar: str, element: str, poscar: str | None = None,
                 emin: float = -10.0, emax: float = 0.0) -> float:
    """One-liner: d-band centre of ``element`` from a DOSCAR."""
    return PDOS(doscar, poscar).dband_center(element, emin=emin, emax=emax)


# --------------------------------------------------------------------------- #
#  Batch scan over a generated tree
# --------------------------------------------------------------------------- #

_LABEL_RE = re.compile(
    r"(?P<host>[A-Za-z]{1,2})(?P<nh>[\d.]+)(?P<dop>[A-Za-z]{1,2})(?P<nd>[\d.]+)")


def parse_label(label: str) -> Dict:
    """Extracts host / dopant / y from a directory label like
    ``Ti2.5V0.5C2O2_s2x2_outer``."""
    head = label.split("_")[0]
    m = _LABEL_RE.match(head)
    out = {"label": label, "pattern": None, "seed": None}
    if m:
        out.update(host=m.group("host"), dopant=m.group("dop"),
                   y=float(m.group("nd")))
    for part in label.split("_")[1:]:
        low = part.lower()
        if part.startswith("s") and "x" in part:
            out["supercell"] = part
        elif low.startswith("seed"):
            out["seed"] = int(low[4:]) if low[4:].isdigit() else None
        elif part.isdigit():
            out["seed"] = int(part)
        elif part in ("outer", "inner", "random", "cluster", "sqs",
                      "random-outer", "random-inner"):
            out["pattern"] = part
    return out


def scan(root: str = "MXenes_alloy", m_elements: Sequence[str] = ("Ti", "V"),
         n: int = 2, emin: float = -10.0, emax: float = 0.0,
         doscar: str = "DOS/DOSCAR", verbose: bool = True
         ) -> List[Dict]:
    """Walks ``root`` and returns one summary dict per completed calculation."""
    rows = []
    target = os.path.basename(doscar)          # os.walk gives BASENAMES only
    for dirpath, _, files in os.walk(root):
        if target not in files:
            continue
        # honour the sub-path too, so --doscar "DOS/DOSCAR" does not also
        # pick up, say, "DOS/PBE0/DOSCAR"
        if os.sep in doscar:
            parent = os.path.basename(os.path.dirname(doscar))
            if os.path.basename(dirpath) != parent:
                continue
        p = os.path.join(dirpath, target)
        try:
            d = PDOS(p)
            row = d.summary(m_elements=m_elements, n=n, emin=emin, emax=emax)
        except Exception as e:                       # noqa: BLE001
            if verbose:
                print(f"  [skip] {dirpath}: {type(e).__name__}: {e}")
            continue
        parts = os.path.relpath(dirpath, root).split(os.sep)
        row.update(parse_label(parts[0] if parts else dirpath))
        row["path"] = dirpath
        rows.append(row)
        if verbose:
            print("  " + "  ".join(f"{k}={v}" for k, v in row.items()
                                   if k in ("label", "eps_d(Ti)", "eps_d(V)")))
    rows.sort(key=lambda r: (r.get("y", 0.0), r.get("pattern") or ""))
    return rows


# --------------------------------------------------------------------------- #
#  Thesis figure
# --------------------------------------------------------------------------- #

def plot_trend(rows: List[Dict], key: str = "eps_d(Ti)",
               extra: Dict[str, Sequence[float]] | None = None,
               out: str = "eps_d_vs_composition.png", dpi: int = 300,
               ylabel: str | None = None, title: str | None = None):
    """eps_d (plus optional experimental series) vs dopant content y.

    ``extra`` maps a series name onto a *list aligned with* ``rows``
    (e.g. ``{"E_ads(H) [eV]": [...], "C_g [F/g]": [...]}``).  Each extra
    series is drawn on a secondary axis after min-max normalisation, so that
    several quantities with different units can be checked for a *common
    trend* -- which is the whole point of the figure.

    Returns ``(path, r_squared)`` of the linear fit of ``key`` vs y.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pts = [(r.get("y"), r.get(key)) for r in rows
           if r.get("y") is not None and r.get(key) is not None]
    pts = [(y, v) for y, v in pts if v == v]         # drop NaN
    if len(pts) < 2:
        raise ValueError(f"Need >= 2 valid points for '{key}', got {len(pts)}")
    pts.sort()
    ys = np.array([p[0] for p in pts], dtype=float)
    vs = np.array([p[1] for p in pts], dtype=float)

    # The extra series are passed by the caller aligned with the *input*
    # `rows`, not with the filtered `pts`.  Restrict them to the same kept
    # rows so plot_trend never indexes with mismatched lengths.
    keep_idx = [i for i, r in enumerate(rows)
                if r.get("y") is not None
                and r.get(key) is not None
                and r.get(key) == r.get(key)]
    extras_aligned: Dict[str, Sequence[float]] = {}
    if extra:
        for name, vals in extra.items():
            extras_aligned[name] = [vals[i] for i in keep_idx]

    coef = np.polyfit(ys, vs, 1)
    fit = np.polyval(coef, ys)
    ss_res = float(((vs - fit) ** 2).sum())
    ss_tot = float(((vs - vs.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.plot(ys, vs, "o", color="#c0392b", ms=6, label=key)
    ax.plot(ys, fit, "-", color="#c0392b", lw=1.2,
            label=f"fit  $R^2$={r2:.3f}")
    ax.set_xlabel(r"V content  $y$  in Ti$_{3-y}$V$_y$C$_2$O$_2$")
    ax.set_ylabel(ylabel or (r"$\varepsilon_d$ (eV, rel. $E_F$)"))

    twin = None
    if extra:
        twin = ax.twinx()
        cmap = plt.cm.Set2(np.linspace(0, 1, max(len(extras_aligned), 3)))
        for k, (name, vals) in enumerate(extras_aligned.items()):
            v = np.array(vals, dtype=float)
            finite = np.isfinite(v)
            # np.ptp(arr) only: ndarray.ptp() was removed in NumPy 2.0
            if finite.sum() < 2 or np.ptp(v[finite]) == 0:
                continue
            vn = (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v))
            twin.plot(ys[finite], vn[finite], "s--", color=cmap[k % len(cmap)],
                      lw=1.0, ms=4, label=name)
        twin.set_ylabel("normalised (0-1)")
        twin.legend(frameon=False, fontsize="x-small", loc="lower right")

    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.legend(frameon=False, fontsize="x-small", loc="upper left")
    if title:
        ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out, r2


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse
    import json

    ap = argparse.ArgumentParser(description="d-band centre analysis.")
    ap.add_argument("doscar", nargs="?", help="a single DOSCAR to analyse")
    ap.add_argument("--poscar", default=None)
    ap.add_argument("--scan", metavar="ROOT", default=None,
                    help="instead: walk a generated tree")
    ap.add_argument("--n", type=int, default=2)
    ap.add_argument("--emin", type=float, default=-10.0)
    ap.add_argument("--emax", type=float, default=0.0)
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("-o", "--out", default="dband_results.json")
    a = ap.parse_args()

    if a.scan:
        rows = scan(a.scan, n=a.n, emin=a.emin, emax=a.emax)
        with open(a.out, "w") as f:
            json.dump(rows, f, indent=2, default=str)
        print(f"\n{len(rows)} calculations -> {a.out}")
        for r in rows:
            print("  " + "  ".join(f"{k}={v}" for k, v in r.items()
                                   if k in ("label", "eps_d(Ti)", "eps_d(V)")))
        if a.plot and len(rows) >= 2:
            for key in ("eps_d(Ti)", "eps_d(V)"):
                if all(key in r for r in rows):
                    p, r2 = plot_trend(
                        rows, key=key,
                        out=key.replace("(", "_").replace(")", "") + "_vs_y.png")
                    print(f"  {key}: R^2 = {r2:.4f}  -> {p}")
    elif a.doscar:
        d = PDOS(a.doscar, a.poscar)
        s = d.summary(n=a.n, emin=a.emin, emax=a.emax)
        print(json.dumps(s, indent=2))
        if a.plot:
            print("plot ->", d.plot(out=a.out if a.out.endswith(".png")
                                    else "pdos.png"))
    else:
        ap.error("give a DOSCAR or --scan ROOT")
