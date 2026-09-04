"""
alloy.py -- Alloyed (multi-element M-site) MXene structure & input generation.

Fork extension for:  "合金化 MXene 电子结构调控与超级电容器储能性能研究"
System:  Ti_(3-y) V_y C2 O2   (n = 2, O-terminated, ABC stacking, H_M hollows)

The upstream `VASP.py` only handles *single*-M MXenes (Mn+1XnT2) in the
primitive 1x1 cell.  This module extends it with:

  1. Multi-element M sublattice (host + dopant) at arbitrary composition y.
  2. Arbitrary in-plane supercells (sx, sy, 1) -- needed to realise
     fractional compositions (e.g. y = 0.5 needs sx*sy even).
  3. Explicit control of *which* M sites the dopant occupies.  For n = 2 the
     three M layers are  M(outer) / M(inner) / M(outer);  the two are NOT
     equivalent, and where the dopant sits is the central question of the
     thesis (V@M1 vs V@M2).
  4. Occupation patterns: ordered (outer / inner), random (seeded,
     reproducible), cluster, and an SQS-lite best-of-N random search that
     maximises pair-disorder.
  5. Corrections to upstream bugs that only appear once M is not a single
     element (POTCAR selection, n inference from the formula).
  6. INCAR generation with DFT+U, vdW (DFT-D3) and the per-species LDAU
     arrays written in POSCAR order (a hard requirement of VASP).

Author of the fork extensions:  Wang Yaoyu (NUIST) -- 2026
Upstream:                       Diego Ontiveros (github.com/diegonti/VASP-MXenes)
"""

from __future__ import annotations

import math
import os
import random
import shutil
from typing import Dict, List, Sequence, Tuple

from VASP import MX, mkdir

# --------------------------------------------------------------------------- #
#  Constants
# --------------------------------------------------------------------------- #

#: Pseudopotential folders that differ from the bare element symbol.
#: (upstream list, kept in sync with VASP.MX.pp)
PP_ALIAS: Dict[str, str] = {
    "Sc": "Sc_sv", "Ti": "Ti_pv", "V": "V_pv", "Cr": "Cr_pv",
    "Y": "Y_sv", "Zr": "Zr_sv", "Nb": "Nb_pv", "Mo": "Mo_pv",
    "Hf": "Hf_pv", "Ta": "Ta_pv", "W": "W_pv",
}

#: Default Hubbard U_eff (eV, Dudarev form) used by the thesis.
DEFAULT_U: Dict[str, float] = {"Ti": 2.5, "V": 2.5, "Nb": 2.0, "Zr": 2.0}

#: Metallic radii (pm), used only for sanity checks / reporting.
METALLIC_RADIUS: Dict[str, float] = {"Ti": 147.0, "V": 134.0, "Nb": 146.0,
                                     "Zr": 160.0, "Mo": 139.0, "Cr": 128.0}

PATTERNS = ("outer", "inner", "random-outer", "random-inner",
            "random", "cluster", "sqs")


# --------------------------------------------------------------------------- #
#  Small helpers
# --------------------------------------------------------------------------- #

def m_site_layers(n: int) -> Tuple[List[int], List[int]]:
    """Outer / inner M-layer indices *within the primitive cell*.

    ``MX.positions()`` emits the M atoms in the order
    ``[M_outer_bottom, M_2, ..., M_(n-1), M_outer_top]``, therefore:

        n = 1  ->  outer = [0, 1], inner = []          (M2X)
        n = 2  ->  outer = [0, 2], inner = [1]         (M3X2)
        n = 3  ->  outer = [0, 3], inner = [1, 2]      (M4X3)
    """
    outer = [0, n]
    inner = list(range(1, n))
    return outer, inner


def suggest_supercell(ys: Sequence[float], max_area: int = 12
                      ) -> Tuple[int, int] | None:
    """Smallest, most square (sx, sy) so that ``sx*sy*y`` is integer for every y.

    A primitive cell contains (n+1) M atoms, so a sx x sy supercell holds
    (n+1)*sx*sy M sites and the dopant count is ``sx*sy*y``.  Fractional y
    therefore requires a commensurate cell:  y = 0.5 needs an even area,
    y = 1.5 likewise.
    """
    best = None
    for area in range(1, max_area + 1):
        for sx in range(1, area + 1):
            if area % sx:
                continue
            sy = area // sx
            if not all(abs(area * y - round(area * y)) < 1e-6 for y in ys):
                continue
            # Smallest area wins FIRST -- DFT cost grows ~N^3, so a 1x3 cell
            # (12 M sites) is ~20x cheaper than a 3x3 one (36 M sites) even
            # though the latter is squarer.  Squareness only breaks ties.
            score = (area, abs(sx - sy))
            if best is None or score < best[0]:
                best = (score, sx, sy)
    return (best[1], best[2]) if best else None


def _fmt(x: float, dec: int = 16) -> str:
    return format(round(float(x), dec), f".{dec}f")


# --------------------------------------------------------------------------- #
#  Occupation patterns
# --------------------------------------------------------------------------- #

def _pair_disorder_score(occ: List[int], neighbours: List[List[int]]) -> float:
    """Lower is more disordered.  Deviation of AA/AB/BB pair counts from random.

    ``occ`` is 0/1 (host/dopant) per site.  Only a cheap heuristic -- for a
    rigorous SQS use ATAT ``mcsqs`` / ICET, but for the 6-27 sites used here
    this is more than sufficient and keeps the fork dependency-free.
    """
    nAA = nAB = nBB = 0
    for i, nbrs in enumerate(neighbours):
        for j in nbrs:
            if j <= i:
                continue
            a, b = occ[i], occ[j]
            if a == 0 and b == 0:
                nAA += 1
            elif a == 1 and b == 1:
                nBB += 1
            else:
                nAB += 1
    n = len(occ)
    x = sum(occ) / n if n else 0.0
    n_pairs = sum(len(v) for v in neighbours) / 2
    # Random-alloy expectation for a symmetric lattice
    expAA = n_pairs * (1 - x) ** 2
    expBB = n_pairs * x ** 2
    expAB = n_pairs * 2 * x * (1 - x)
    return ((nAA - expAA) ** 2 + (nBB - expBB) ** 2 + (nAB - expAB) ** 2) ** 0.5


def _neighbour_list(sites: List[List[float]], cutoff: float = 1.05,
                    lattice: Tuple[float, float] = (1.0, 1.0)) -> List[List[int]]:
    """Minimum-image in-plane neighbour list for fractional (a, b) coordinates.

    ``cutoff`` is in units where the in-plane M-M distance of the primitive
    hexagonal cell is ~0.577 (1/sqrt(3)).  The default therefore grabs the
    6 first in-plane neighbours of each M site.
    """
    nbrs = [[] for _ in sites]
    n = len(sites)
    for i in range(n):
        ax, ay = sites[i][0], sites[i][1]
        for j in range(n):
            if i == j:
                continue
            dx = sites[j][0] - ax
            dy = sites[j][1] - ay
            dx -= round(dx)
            dy -= round(dy)
            d = math.hypot(dx * lattice[0], dy * lattice[1])
            if 1e-6 < d <= cutoff:
                nbrs[i].append(j)
    return nbrs


def _round_robin(groups: List[List[int]], n_dopant: int,
                 rng: random.Random | None = None) -> List[int]:
    """Picks ``n_dopant`` sites by taking one from each layer group in turn.

    This is what keeps the two *outer* surfaces of the slab equivalent.  The
    naive "fill the first N outer sites" puts every dopant on the bottom
    surface only, which creates a spurious slab dipole and poisons both the
    total energy and the work function.
    """
    if rng is not None:
        groups = [list(g) for g in groups]
        for g in groups:
            rng.shuffle(g)
    picked: List[int] = []
    cursor = [0] * len(groups)
    while len(picked) < n_dopant:
        progressed = False
        for gi, g in enumerate(groups):
            if cursor[gi] < len(g) and len(picked) < n_dopant:
                picked.append(g[cursor[gi]])
                cursor[gi] += 1
                progressed = True
        if not progressed:                       # every group exhausted
            raise ValueError(f"cannot place {n_dopant} dopants in "
                             f"{sum(len(g) for g in groups)} available sites")
    return picked


def build_occupation(n_sites: int,
                     layer_of_site: List[str],
                     n_dopant: int,
                     pattern: str = "outer",
                     seed: int = 42,
                     sites_frac: List[List[float]] | None = None,
                     groups: Dict[int, List[int]] | None = None,
                     n_random_tries: int = 400) -> List[int]:
    """Returns a 0/1 list (0 = host, 1 = dopant) of length ``n_sites``.

    Parameters
    ----------
    n_sites : total number of M sites in the supercell.
    layer_of_site : ``"outer"`` / ``"inner"`` tag per site.
    n_dopant : how many sites the dopant occupies.
    pattern : one of ``PATTERNS``.
    seed : RNG seed -- every pattern is deterministic, so a case is
        reproducible from its config alone (essential for a thesis).
    groups : ``{primitive_M_layer_index: [site indices]}``.  When supplied,
        the ordered patterns distribute the dopants round-robin over the
        layer groups so the slab stays symmetric.
    """
    if n_dopant > n_sites:
        raise ValueError(f"n_dopant={n_dopant} > n_sites={n_sites}")
    if pattern not in PATTERNS:
        raise ValueError(f"pattern must be one of {PATTERNS}, got {pattern!r}")

    occ = [0] * n_sites
    if n_dopant == 0:
        return occ
    if n_dopant == n_sites:
        return [1] * n_sites

    outer_idx = [i for i, l in enumerate(layer_of_site) if l == "outer"]
    inner_idx = [i for i, l in enumerate(layer_of_site) if l == "inner"]

    # Layer groups restricted to the requested sublattice, kept in a stable
    # order (bottom-outer, inner..., top-outer).
    def layer_groups(which: str) -> List[List[int]]:
        if groups is None:
            return [[i for i in (outer_idx if which == "outer" else inner_idx)]]
        return [sorted(g) for k, g in sorted(groups.items())
                if g and layer_of_site[g[0]] == which]

    # ---------------- ordered (deterministic, round-robin over layers) ----
    if pattern in ("outer", "inner"):
        gs = layer_groups(pattern)
        have = sum(len(g) for g in gs)
        if n_dopant > have:
            raise ValueError(f"pattern {pattern!r} cannot host {n_dopant} "
                             f"dopants (only {have} {pattern} sites)")
        for i in _round_robin(gs, n_dopant):
            occ[i] = 1
        return occ

    # ---------------- cluster (deliberately bad, used as a control) -------
    if pattern == "cluster":
        start = n_sites // 2 - n_dopant // 2
        for i in range(start, start + n_dopant):
            occ[i % n_sites] = 1
        return occ

    # ---------------- random family --------------------------------------
    if pattern == "random-outer":
        gs = layer_groups("outer")
    elif pattern == "random-inner":
        gs = layer_groups("inner")
    else:                                            # "random" / "sqs"
        gs = layer_groups("outer") + layer_groups("inner")

    if n_dopant > sum(len(g) for g in gs):
        raise ValueError(f"pattern {pattern!r} cannot host {n_dopant} dopants "
                         f"({sum(len(g) for g in gs)} sites available)")

    rng = random.Random(seed)

    if pattern != "sqs" or sites_frac is None:
        for i in _round_robin(gs, n_dopant, rng=rng):
            occ[i] = 1
        return occ

    # ---- SQS-lite: best of N balanced-random draws by pair-disorder ------
    # Candidates are generated with the *same* round-robin balancing, so the
    # search only optimises the in-plane arrangement, not the layer balance.
    nbrs = _neighbour_list(sites_frac)
    best_occ, best_score = None, None
    for _ in range(n_random_tries):
        cand = [0] * n_sites
        for i in _round_robin(gs, n_dopant, rng=rng):
            cand[i] = 1
        s = _pair_disorder_score(cand, nbrs)
        if best_score is None or s < best_score - 1e-12:
            best_occ, best_score = cand, s
    return best_occ


# --------------------------------------------------------------------------- #
#  The alloyed MXene
# --------------------------------------------------------------------------- #

class AlloyMX(MX):
    """Alloyed MXene  (Host_(n+1-y) Dopant_y) X_n T_2  in a sx x sy supercell.

    Parameters
    ----------
    host, dopant : element symbols of the M sublattice (e.g. "Ti", "V").
    y : dopant count *per M(n+1) formula unit* (0 .. n+1).
        y = 0.5 means Ti2.5V0.5C2 for n = 2.
    n : MXene width (2 -> M3X2, the thesis system).
    X : "C" or "N".
    T : termination, "" (pristine), "O2", "F2", "O2H2", ...
    supercell : (sx, sy, 1).
    pattern : dopant site-occupation pattern (see :data:`PATTERNS`).
    seed : RNG seed for the random patterns.
    """

    # ---------------------------------------------------------------- setup --
    def __init__(self, host: str = "Ti", dopant: str = "V", y: float = 0.0,
                 n: int = 2, X: str = "C", T: str = "O2",
                 supercell: Tuple[int, int, int] = (2, 2, 1),
                 pattern: str = "outer", seed: int = 42,
                 stacking: str = "ABC", hollows: str = "HM",
                 do: float = 2.5, v: float = 20.0, ab: float = 3.0,
                 lp: float = 1.0, root: str = "MXenes_alloy"):

        self.host, self.dopant = host, dopant
        self.y = float(y)
        self.n = int(n)
        self.Xel = X
        self.T = T
        self.pattern = pattern
        self.seed = int(seed)
        self.root = root

        sx, sy, sz = supercell
        self.sx, self.sy, self.sz = int(sx), int(sy), int(sz)
        self.area = self.sx * self.sy

        # ---------------- composition bookkeeping -------------------------- #
        nM_prim = self.n + 1                       # 3 for M3X2
        self.nM_prim = nM_prim
        self.nM_sites = nM_prim * self.area
        self.nX_sites = self.n * self.area
        self.nT_sites = 0 if not T else 2 * self.area

        n_dop_float = self.area * self.y
        self.n_dopant = int(round(n_dop_float))
        if abs(n_dop_float - self.n_dopant) > 1e-6:
            raise ValueError(
                f"Composition y={y} is not commensurate with a "
                f"{self.sx}x{self.sy} supercell (needs {n_dop_float} dopants). "
                f"Use suggest_supercell() to pick a valid cell.")
        self.n_host = self.nM_sites - self.n_dopant

        # ---------------- geometry ---------------------------------------- #
        self.stacking = stacking
        self.hollows = hollows
        self.do = do
        self.v = v
        self.ab = ab
        self.lp = lp
        self.shift = 1 if T else 0
        self.terminal = bool(T)
        self.pristine = f"{host}{nM_prim}{X}{self.n}"

        # Upstream MX writes/reads these; keep them consistent.
        self.mxName = (f"{host}{nM_prim - self.y:g}{dopant}{self.y:g}"
                       f"{X}{self.n}{T}")
        self.name = self.mxName
        self.label = (f"{self.mxName}_s{self.sx}x{self.sy}_{pattern}"
                      + (f"_seed{seed}" if "random" in pattern or pattern == "sqs"
                         else ""))
        self.atoms = [a for a in (host, dopant, X, self._t_atoms()) if a]
        self.index = [c for c in (self.n_host, self.n_dopant,
                                  self.nX_sites, self.nT_sites) if c]

        # directory tree
        self.pdir = f"./{root}/{self.label}/{self.stacking}/{self.hollows}/"
        calc_dirs = ["opt", "DOS", "BS", "WF", "bader"]
        self.dirs = [self.pdir + d + "/" for d in calc_dirs]
        self.aimsBS_dirs = []

        self._U: Dict[str, float] = dict(DEFAULT_U)

    def _t_atoms(self) -> str:
        """Element symbol of the termination ("" / "O" / "F" / "OH"->"O"+"H")."""
        if not self.T:
            return ""
        body = "".join(ch for ch in self.T if not ch.isdigit())
        return body if body in ("O", "F", "Cl", "S", "Br", "I") else ""

    # ------------------------------------------------------------ geometry --
    def primitive_positions(self):
        """(M, X, T) fractional positions of the *primitive* cell.

        Delegates to upstream ``MX.positions()`` (which is well tested for
        ABC / ABA stackings and all hollow sites) via a throw-away MX object.
        """
        prim_name = f"{self.host}{self.nM_prim}{self.Xel}{self.n}{self.T}"
        base = MX(prim_name, self.stacking, self.hollows)
        base.n, base.do, base.v, base.ab, base.lp = self.n, self.do, self.v, self.ab, self.lp
        base.shift = self.shift
        M, X, T = base.positions()
        # Upstream returns pre-formatted *strings*; the supercell arithmetic
        # needs real numbers, so cast back here.
        to_f = lambda block: [[float(c) for c in p] for p in block]
        return to_f(M), to_f(X), to_f(T)

    def build(self) -> Tuple[List[Tuple[str, List[float]]], List[str], List[int]]:
        """Builds the supercell.

        Returns
        -------
        sites : list of ``(element, [a, b, c])`` for **every** atom, in the
            canonical order  M_outer(bottom) -> X -> M_inner -> ... -> T,
            replicated over the supercell.
        species : POSCAR species line, e.g. ``["Ti", "V", "C", "O"]``.
        counts : POSCAR counts line, matching ``species``.
        """
        M, X, T = self.primitive_positions()

        # --- replicate the primitive cell over sx x sy --------------------
        # M_all is built "layer-major": [M_k=0 for every cell, M_k=1 for every
        # cell, ...]  X_all / T_all the same way.  `m_index[(k, cell)]` maps
        # (primitive M layer k, cell replica c) onto its position in M_all.
        outer_i, inner_i = m_site_layers(self.n)

        M_all: List[List[float]] = []
        layer_of_site: List[str] = []
        m_index: Dict[Tuple[int, int], int] = {}
        layer_groups: Dict[int, List[int]] = {}
        for k in range(len(M)):
            layer = "outer" if k in outer_i else "inner"
            layer_groups[k] = []
            for c in range(self.area):
                da = (c % self.sx) / self.sx
                db = (c // self.sx) / self.sy
                M_all.append([M[k][0] / self.sx + da,
                              M[k][1] / self.sy + db,
                              M[k][2]])
                layer_of_site.append(layer)
                m_index[(k, c)] = len(M_all) - 1
                layer_groups[k].append(len(M_all) - 1)

        def _replicate(prim):
            out = []
            for p in prim:
                for c in range(self.area):
                    da = (c % self.sx) / self.sx
                    db = (c // self.sx) / self.sy
                    out.append([p[0] / self.sx + da,
                                p[1] / self.sy + db, p[2]])
            return out

        X_all = _replicate(X)
        T_all = _replicate(T)

        # --- decide which M sites carry the dopant -------------------------
        if self.n_dopant == 0:
            occ = [0] * len(M_all)
        elif self.n_dopant == len(M_all):
            occ = [1] * len(M_all)
        else:
            occ = build_occupation(len(M_all), layer_of_site,
                                   self.n_dopant, self.pattern, self.seed,
                                   sites_frac=M_all, groups=layer_groups)

        self.occ = occ
        self.layer_of_site = layer_of_site
        self.layer_groups = layer_groups
        self.m_index = m_index

        # --- reassemble the slab: M0 X0 M1 X1 ... Mn  (then T) -------------
        sites: List[Tuple[str, List[float]]] = []
        for c in range(self.area):
            for k in range(self.nM_prim):
                s = m_index[(k, c)]
                el = self.dopant if occ[s] else self.host
                sites.append((el, M_all[s]))
                if k < self.n:
                    sites.append((self.Xel, X_all[k * self.area + c]))
        for p in T_all:
            sites.append((self._t_atoms(), p))

        # --- group by species (VASP requires contiguous blocks) ------------
        # Canonical order  host -> dopant -> X -> T  (species with zero atoms
        # are dropped), so that the POTCAR / INCAR ordering stays identical
        # across every composition and pattern.  Deriving the order from first
        # appearance would silently shuffle it whenever y hits 0 or n+1.
        canonical = [self.host, self.dopant, self.Xel, self._t_atoms()]
        order = [el for el in dict.fromkeys(canonical)
                 if el and any(e == el for e, _ in sites)]
        grouped = [(el, p) for el in order for (e, p) in sites if e == el]
        counts = [sum(1 for e, _ in sites if e == el) for el in order]

        self.species = order
        self.counts = counts
        self.nAtoms = sum(counts)
        return grouped, order, counts

    # ------------------------------------------------------------- POSCAR --
    def POSCAR(self, sd_flags: str = "T T T", comment: str | None = None):
        """Writes POSCAR into every ``opt/`` directory."""
        grouped, species, counts = self.build()
        a = self.ab * self.sx
        b = self.ab * self.sy
        c = self.n * self.do + 2 * self.shift + self.v

        f1 = ".10f"
        lat_a = [format(a, f1), format(0, f1), format(0, f1)]
        lat_b = [format(round(-b * math.sin(math.pi / 6), 10), f1),
                 format(round(b * math.cos(math.pi / 6), 10), f1),
                 format(0, f1)]
        lat_c = [format(0, f1), format(0, f1), format(c, f1)]

        body = ""
        for el, pos in grouped:
            body += (f"  {_fmt(pos[0])}  {_fmt(pos[1])}  {_fmt(pos[2])}"
                     f"  {sd_flags}\n")

        header = comment or (f"{self.mxName} {self.sx}x{self.sy} "
                             f"{self.pattern} {self.stacking}-{self.hollows}")

        for d in self.dirs:
            if not d.endswith("opt/"):
                continue
            with open(d + "POSCAR", "w") as f:
                f.write(f"{header}\n"
                        f"   {format(self.lp, '.14f')}\n"
                        f"     {lat_a[0]}   {lat_a[1]}   {lat_a[2]}\n"
                        f"    {lat_b[0]}   {lat_b[1]}   {lat_b[2]}\n"
                        f"     {lat_c[0]}   {lat_c[1]}  {lat_c[2]}\n"
                        f"  {'  '.join(species)}\n"
                        f"  {'  '.join(str(i) for i in counts)}\n"
                        f"Selective Dynamics\nDirect\n{body}")

    # ------------------------------------------------------------- POTCAR --
    def POTCAR(self):
        """Concatenates POTCARs **in POSCAR species order**.

        Upstream ``MX.POTCAR`` only special-cases the *first* atom
        (``if self.atoms[0] in p and first``), so as soon as a second metal is
        present the dopant silently falls back to the *standard* PAW potential
        (V instead of V_pv).  That is exactly the alloy case, so it is fixed
        here: every species is mapped through :data:`PP_ALIAS`.
        """
        for d in [self.pdir] + self.dirs:
            with open(d + "POTCAR", "w") as out:
                for el in self.species:
                    folder = PP_ALIAS.get(el, el)
                    path = f"./PP/{folder}/POTCAR"
                    if not os.path.exists(path):
                        raise FileNotFoundError(
                            f"POTCAR for {el} not found at {path}. "
                            f"Add the pseudopotential or edit PP_ALIAS.")
                    with open(path, "r") as src:
                        out.write(src.read())

    # ------------------------------------------------------------ KPOINTS --
    def KPOINTS(self, p1: int = 7, p2: int = 7, p3: int = 1,
                mesh: Tuple[int, int, int] | None = None):
        """Gamma-centred mesh, automatically reduced for supercells.

        For an sx x sy supercell the primitive mesh (p1, p2) is divided by
        (sx, sy) so the *reciprocal-space density* stays constant -- running a
        2x2 cell with 7x7x1 is wasteful and slow.
        """
        if mesh is not None:
            k1, k2, k3 = mesh
        else:
            k1 = max(1, int(round(p1 / self.sx)))
            k2 = max(1, int(round(p2 / self.sy)))
            k3 = p3
        self.kmesh = (k1, k2, k3)

        txt = (f"K-Points\n0\nGamma\n {k1}  {k2}  {k3}\n 0  0  0\n")
        for d in self.dirs:
            with open(d + "KPOINTS", "w") as f:
                f.write(txt)

    # -------------------------------------------------------------- INCAR --
    def _ldau_block(self) -> str:
        """LDAU / LDAUU / LDAUJ / LDAUL arrays in POSCAR species order."""
        if not self._U:
            return ""
        lda_u, lda_j, lda_l = [], [], []
        for el in self.species:
            u = self._U.get(el, 0.0)
            lda_u.append(f"{u:.2f}")
            lda_j.append("0.00")
            lda_l.append("2" if el in DEFAULT_U else "-1")
        return (f"LDAU   = .TRUE.\n"
                f"LDAUTYPE = 2\n"
                f"LDAUPRINT = 2\n"
                f"LDAUU  = {'  '.join(lda_u)}\n"
                f"LDAUJ  = {'  '.join(lda_j)}\n"
                f"LDAUL  = {'  '.join(lda_l)}\n")

    def set_U(self, **kwargs):
        """Override Hubbard U per element, e.g. ``mx.set_U(Ti=2.5, V=2.5)``."""
        self._U.update({k: float(v) for k, v in kwargs.items()})
        return self

    def INCAR(self, encut: float = 550, vdw: str = "d3", spin: bool = True,
              ediff: str = "1E-06", ediffg: str = "-0.01",
              ismear: float | None = None, sigma: float = 0.05,
              lorbit: int = 11, nedos: int = 5000,
              dipol: bool = True):
        """Writes a per-calculation INCAR.

        Notes on the choices (these follow the experimental plan of the thesis):
          * ``encut = 550`` eV  -- converged for Ti/V/C/O PAW sets.
          * ``vdw = "d3"``      -- Grimme D3, mandatory for layered MXenes.
          * ``LDAU``            -- Dudarev U_eff, 2.5 eV for both Ti and V.
          * ``LDIPOL/IDIPOL=3`` -- slab dipole correction, required for a
            meaningful work function (LOCPOT) and for any asymmetric
            adsorption case.
          * ``LORBIT = 11``     -- needed for the lm-decomposed DOSCAR that
            :mod:`dband` integrates to get the d-band centre.
        """
        per_calc = {
            "opt":    dict(IBRION=2, NSW=200, ISIF=2, ISMEAR=0, LWAVE=".FALSE."),
            "DOS":    dict(IBRION=-1, NSW=0, ISMEAR=-5, LORBIT=lorbit,
                           EMIN=-15, EMAX=15, NEDOS=nedos, LWAVE=".FALSE."),
            "BS":     dict(IBRION=-1, NSW=0, ISMEAR=0, LORBIT=lorbit,
                           LWAVE=".TRUE.", ICHARG=11),
            "WF":     dict(IBRION=-1, NSW=0, ISMEAR=-5, LORBIT=lorbit,
                           LVTOT=".TRUE.", LVHAR=".TRUE.", NEDOS=nedos,
                           NGXF=120, NGYF=120, NGZF=700),
            "bader":  dict(IBRION=-1, NSW=0, ISMEAR=-5, LORBIT=lorbit,
                           LAECHG=".TRUE.", LVTOT=".TRUE.", NEDOS=nedos),
        }
        vdw_map = {"d3": 11, "d2": 10, "d3bj": 12, "optb86b": 3, "none": None}

        for d in self.dirs:
            key = d.rstrip("/").split("/")[-1]
            extra = per_calc.get(key, {})
            sm = extra.pop("ISMEAR") if "ISMEAR" in extra else ismear
            if sm is None:
                sm = 0

            lines = [
                f"SYSTEM = {self.label}_{key}",
                "",
                "# ---------- electronic ----------",
                f"ENCUT  = {encut}",
                f"EDIFF  = {ediff}",
                f"EDIFFG = {ediffg}",
                f"ISPIN  = {2 if spin else 1}",
                f"ISMEAR = {sm}",
                f"SIGMA  = {sigma}",
                "LREAL  = AUTO",
                "ALGO   = Normal",
                "NELMIN = 4",
                "NELM   = 200",
                "PREC   = Accurate",
                "GGA    = PE",
                "",
                "# ---------- Hubbard / vdW ----------",
            ]
            if self._U:
                lines += ["LDAU   = .TRUE.", "LASPH  = .TRUE.",
                          "LDAUTYPE = 2", "LDAUPRINT = 2",
                          "LDAUU  = " + "  ".join(
                              f"{self._U.get(e, 0.0):.2f}" for e in self.species),
                          "LDAUJ  = " + "  ".join("0.00" for _ in self.species),
                          "LDAUL  = " + "  ".join(
                              "2" if e in DEFAULT_U else "-1" for e in self.species)]
            if vdw_map.get(vdw):
                lines.append(f"IVDW   = {vdw_map[vdw]}")
            lines += ["",
                      "# ---------- ionic / IO ----------",
                      "ISTART = 0", "INIWAV = 1", "NWRITE = 2",
                      "LCHARG = .FALSE."]
            if spin:
                # MAGMOM is per ATOM (NIONS entries), not per species -- the
                # VASP "n*value" shorthand keeps it readable for supercells.
                mag = [f"{c}*{'1.0' if e in DEFAULT_U else '0.6'}"
                       for e, c in zip(self.species, self.counts)]
                lines.append("MAGMOM = " + "  ".join(mag))
            lines += ["", "# ---------- calculation specific ----------"]
            for k, v in extra.items():
                lines.append(f"{k:7s} = {v}")
            if dipol:
                lines += ["LDIPOL = .TRUE.", "IDIPOL = 3", "DIPOL  = 0.5 0.5 0.5"]
            lines += ["", "# " + "-" * 40,
                      f"# {self.mxName}  supercell {self.sx}x{self.sy}x{self.sz}"
                      f"  pattern={self.pattern}  seed={self.seed}",
                      f"# species order: {' '.join(self.species)}"
                      f"  counts: {' '.join(map(str, self.counts))}", ""]

            with open(d + "INCAR", "w") as f:
                f.write("\n".join(lines) + "\n")

    # -------------------------------------------------------------- driver --
    def createDirs(self):
        mkdir(f"./{self.root}/",
              f"./{self.root}/{self.label}/",
              f"./{self.root}/{self.label}/{self.stacking}/",
              f"./{self.root}/{self.label}/{self.stacking}/{self.hollows}/")
        mkdir(*self.dirs)

    def write_inputs(self, script: str | None = None, **incar_kw):
        """Creates the tree and writes POSCAR / POTCAR / KPOINTS / INCAR.

        ``POSCAR()`` also populates ``self.species`` / ``self.counts``, which
        POTCAR and INCAR depend on -- so the order of calls matters.
        """
        self.createDirs()
        self.POSCAR()          # -> self.species, self.counts
        self.POTCAR()
        self.KPOINTS()
        self.INCAR(**incar_kw)
        if script and os.path.exists(script):
            for d in self.dirs:
                shutil.copy(script, d + "script")
        return self

    # ---------------------------------------------------------- reporting --
    def composition(self) -> Dict[str, float]:
        """Composition normalised to one M(n+1)XnT2 formula unit."""
        f = self.area
        out = {self.host: round(self.n_host / f, 4),
               self.dopant: round(self.n_dopant / f, 4),
               self.Xel: float(self.n),
               self._t_atoms(): 2.0}
        return {k: v for k, v in out.items() if k and v}

    def site_summary(self) -> Dict[str, int]:
        """How many dopant atoms sit on outer vs inner M layers."""
        if not hasattr(self, "occ"):
            self.build()
        o = sum(1 for i, x in enumerate(self.occ)
                if x == 1 and self.layer_of_site[i] == "outer")
        i_ = sum(self.occ) - o
        return {"dopant_outer": o, "dopant_inner": i_,
                "total_dopant": sum(self.occ),
                "n_outer_sites": self.layer_of_site.count("outer"),
                "n_inner_sites": self.layer_of_site.count("inner")}

    def __repr__(self):
        return (f"AlloyMX({self.mxName}, {self.sx}x{self.sy}, "
                f"pattern={self.pattern}, atoms={self.nAtoms})")


# --------------------------------------------------------------------------- #
#  Batch generation
# --------------------------------------------------------------------------- #

def generate_batch(ys: Sequence[float], patterns: Sequence[str] = ("outer", "inner"),
                   host: str = "Ti", dopant: str = "V", n: int = 2,
                   X: str = "C", T: str = "O2", supercell=(2, 2, 1),
                   seed: int = 42, stacking: str = "ABC", hollows: str = "HM",
                   v: float = 20.0, root: str = "MXenes_alloy",
                   script: str | None = None, verbose: bool = True,
                   **incar_kw) -> List[AlloyMX]:
    """Generates the full composition x pattern matrix the thesis needs."""
    made = []
    for y in ys:
        for pat in patterns:
            mx = AlloyMX(host=host, dopant=dopant, y=y, n=n, X=X, T=T,
                         supercell=supercell, pattern=pat, seed=seed,
                         stacking=stacking, hollows=hollows, v=v, root=root)
            try:
                mx.write_inputs(script=script, **incar_kw)
            except ValueError as e:
                if verbose:
                    print(f"  [skip] {mx.mxName} {pat}: {e}")
                continue
            made.append(mx)
            if verbose:
                s = mx.site_summary()
                print(f"  {mx.label:<46s} atoms={mx.nAtoms:<4d} "
                      f"V(outer/inner)={s['dopant_outer']}/{s['dopant_inner']} "
                      f"k={mx.kmesh}")
    return made


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Generate VASP inputs for alloyed MXenes "
                    "Ti_(3-y)V_yC2O2 and friends.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--host", default="Ti")
    p.add_argument("--dopant", default="V")
    p.add_argument("-y", "--ys", nargs="+", type=float,
                   default=[0.0, 0.5, 1.0, 1.5],
                   help="dopant count per M(n+1) formula unit")
    p.add_argument("-n", type=int, default=2, help="MXene width (2 -> M3X2)")
    p.add_argument("-X", default="C")
    p.add_argument("-T", default="O2")
    p.add_argument("-s", "--supercell", nargs=2, type=int, default=[2, 2])
    p.add_argument("-p", "--patterns", nargs="+", default=["outer", "inner", "sqs"])
    p.add_argument("--stacking", default="ABC")
    p.add_argument("--hollows", default="HM")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("-v", "--vacuum", type=float, default=20.0)
    p.add_argument("--root", default="MXenes_alloy")
    p.add_argument("--script", default=None)
    p.add_argument("--encut", type=float, default=550)
    p.add_argument("--vdw", default="d3")
    a = p.parse_args()

    print(f"Suggesting supercell for y={a.ys}: "
          f"{suggest_supercell(a.ys)}")
    print(f"Generating {len(a.ys)} compositions x {len(a.patterns)} patterns\n")

    made = generate_batch(a.ys, a.patterns, host=a.host, dopant=a.dopant,
                          n=a.n, X=a.X, T=a.T,
                          supercell=(a.supercell[0], a.supercell[1], 1),
                          seed=a.seed, stacking=a.stacking, hollows=a.hollows,
                          v=a.vacuum, root=a.root, script=a.script,
                          encut=a.encut, vdw=a.vdw)
    print(f"\nDone: {len(made)} cases written under ./{a.root}/")
