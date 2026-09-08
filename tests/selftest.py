"""
selftest.py -- end-to-end verification of the alloy-MXene extensions.

Everything here runs **without VASP**: structures are generated from scratch
and the DOSCAR analysed is synthetic, with an analytically known d-band
centre.  The point is to catch parser/geometry regressions before they cost
a week of machine time.

    python tests/selftest.py            # run and clean up
    python tests/selftest.py -v         # keep the work directory
    python -m pytest tests/selftest.py  # also works under pytest

Note: POTCAR lookup is relative to the CWD (``./PP/...``), so the tests run
inside the repository root and generate into ``_selftest_work/``.
"""

from __future__ import annotations

import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
os.chdir(ROOT)                      # ./PP must resolve for POTCAR()

from alloy import AlloyMX, suggest_supercell, m_site_layers          # noqa: E402
from dband import PDOS, read_poscar, parse_label                     # noqa: E402

WORK = "_selftest_work"
PASSED: list = []
FAILED: list = []


def check(name: str, cond: bool, detail: str = ""):
    (PASSED if cond else FAILED).append(name)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" +
          (f"   {detail}" if detail else ""))


def approx(a, b, tol=1e-3):
    return abs(float(a) - float(b)) <= tol


def parse_incar(path: str) -> dict:
    """Parses an INCAR into ``{KEY: value}``, ignoring comments/formatting.

    Testing on the parsed dict rather than on raw strings means the checks
    stay valid if the INCAR writer's whitespace ever changes.
    """
    out = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        k, v = k.strip(), v.strip()
        if k and v:
            out[k.upper()] = v
    return out


def _slab_poscar(pattern: str, y: float = 1.0) -> str:
    """Generates (once) an alloy POSCAR and returns its path."""
    mx = AlloyMX(host="Ti", dopant="V", y=y, n=2, T="O2",
                 supercell=(2, 2, 1), pattern=pattern, root=WORK)
    p = os.path.join(WORK, mx.label, "ABC", "HM", "opt", "POSCAR")
    if not os.path.exists(p):
        mx.write_inputs(encut=550, vdw="d3")
    return p


# --------------------------------------------------------------------------- #
#  1. alloy.py -- composition, supercell, site balance
# --------------------------------------------------------------------------- #

def test_alloy_basics():
    print("\n[1] alloy.py -- structure generation")

    # Smallest area wins: a 1x2 cell (area 2, 6 M sites) is commensurate with
    # y = 0.5/1.0/1.5 and is ~2x cheaper than 2x2.  Squareness only breaks
    # ties, never overrides cost.
    check("suggest_supercell([0,0.5,1,1.5]) == (1,2)",
          suggest_supercell([0.0, 0.5, 1.0, 1.5]) == (1, 2),
          str(suggest_supercell([0.0, 0.5, 1.0, 1.5])))
    check("suggest_supercell([0,1/3]) is (1,3) or (3,1)",
          suggest_supercell([0.0, 1 / 3]) in ((1, 3), (3, 1)),
          str(suggest_supercell([0.0, 1 / 3])))
    check("m_site_layers(2) -> outer=[0,2] inner=[1]",
          m_site_layers(2) == ([0, 2], [1]))

    mx = AlloyMX(host="Ti", dopant="V", y=1.0, n=2, T="O2",
                 supercell=(2, 2, 1), pattern="outer")
    _, species, counts = mx.build()

    check("M3C2O2 in a 2x2 cell has 28 atoms", mx.nAtoms == 28, f"got {mx.nAtoms}")
    check("atom counts Ti8 V4 C8 O8", counts == [8, 4, 8, 8], f"got {counts}")
    check("species order is host,dopant,X,T", species == ["Ti", "V", "C", "O"],
          f"got {species}")
    check("composition normalises to Ti2V1C2O2",
          mx.composition() == {"Ti": 2.0, "V": 1.0, "C": 2.0, "O": 2.0},
          f"got {mx.composition()}")

    try:
        AlloyMX(y=0.5, supercell=(3, 3, 1), pattern="outer")
        check("non-commensurate supercell rejected", False)
    except ValueError:
        check("non-commensurate supercell rejected", True)

    mx0 = AlloyMX(y=0.0, pattern="outer"); mx0.build()
    check("y=0 drops the dopant species", mx0.species == ["Ti", "C", "O"],
          f"got {mx0.species}")
    mx3 = AlloyMX(y=3.0, pattern="outer"); mx3.build()
    check("y=3 (hypothetical V3C2O2 end member) works",
          mx3.species == ["V", "C", "O"], f"got {mx3.species}")


def test_alloy_symmetry():
    print("\n[2] alloy.py -- surface symmetry (no spurious slab dipole)")

    for y in (0.5, 1.0, 1.5):
        mx = AlloyMX(y=y, pattern="outer"); mx.build()
        groups = [g for k, g in sorted(mx.layer_groups.items())
                  if mx.layer_of_site[g[0]] == "outer"]
        per_layer = [sum(mx.occ[i] for i in g) for g in groups]
        ok = len(per_layer) == 2 and abs(per_layer[0] - per_layer[1]) <= 1
        check(f"y={y}: V balanced across the two outer surfaces", ok,
              f"bottom={per_layer[0]} top={per_layer[1]}")

    mx = AlloyMX(y=1.0, pattern="inner"); mx.build()
    s = mx.site_summary()
    check("pattern=inner puts every V on the inner layer",
          s["dopant_outer"] == 0 and s["dopant_inner"] == 4, str(s))

    a = AlloyMX(y=1.0, pattern="sqs", seed=7); a.build()
    b = AlloyMX(y=1.0, pattern="sqs", seed=7); b.build()
    check("sqs is reproducible for a fixed seed", a.occ == b.occ)


def test_alloy_inputs():
    print("\n[3] alloy.py -- written VASP inputs")

    mx = AlloyMX(host="Ti", dopant="V", y=1.0, n=2, T="O2",
                 supercell=(2, 2, 1), pattern="outer", root=WORK)
    mx.write_inputs(encut=550, vdw="d3")
    d = os.path.join(WORK, mx.label, "ABC", "HM")

    pos = open(os.path.join(d, "opt", "POSCAR")).read().splitlines()
    check("POSCAR species line", pos[5].split() == ["Ti", "V", "C", "O"], pos[5])
    check("POSCAR counts line",
          [int(x) for x in pos[6].split()] == [8, 4, 8, 8], pos[6])
    check("POSCAR has 28 coordinate lines", len(pos) == 9 + 28, str(len(pos)))

    pot = open(os.path.join(d, "opt", "POTCAR")).read()
    check("POTCAR uses Ti_pv (not bare Ti)", "PAW_PBE Ti_pv" in pot)
    check("POTCAR uses V_pv (upstream fell back to bare V)", "PAW_PBE V_pv" in pot)
    check("POTCAR order matches POSCAR (Ti before V)",
          pot.find("Ti_pv") < pot.find("V_pv"))

    inc = parse_incar(os.path.join(d, "opt", "INCAR"))
    check("INCAR ENCUT = 550", inc.get("ENCUT") == "550")
    check("INCAR IVDW = 11 (DFT-D3)", inc.get("IVDW") == "11")
    check("INCAR LDAUU with U(Ti)=U(V)=2.5",
          inc.get("LDAUU") == "2.50  2.50  0.00  0.00", str(inc.get("LDAUU")))
    check("INCAR LDAUL = 2 for both metals",
          inc.get("LDAUL") == "2  2  -1  -1", str(inc.get("LDAUL")))
    check("INCAR slab dipole correction", inc.get("IDIPOL") == "3")
    check("INCAR MAGMOM expanded per atom",
          inc.get("MAGMOM") == "8*1.0  4*1.0  8*0.6  8*0.6", str(inc.get("MAGMOM")))
    check("DOS INCAR sets LORBIT=11 (required by dband.py)",
          parse_incar(os.path.join(d, "DOS", "INCAR")).get("LORBIT") == "11")
    check("KPOINTS reduced for the 2x2 cell (4 4 1)",
          " 4  4  1" in open(os.path.join(d, "DOS", "KPOINTS")).read())


# --------------------------------------------------------------------------- #
#  2. dband.py -- synthetic DOSCAR with a known answer
# --------------------------------------------------------------------------- #

def make_synthetic(case: str, centres: dict, nedos: int = 2001):
    """Writes DOSCAR + POSCAR whose d-band centres are exactly ``centres``."""
    mx = AlloyMX(host="Ti", dopant="V", y=1.0, n=2, T="O2",
                 supercell=(2, 2, 1), pattern="outer", root=WORK)
    mx.write_inputs()
    src = os.path.join(WORK, mx.label, "ABC", "HM", "opt", "POSCAR")

    E = np.linspace(-10, 6, nedos)
    EF = 3.0
    Eabs = E + EF
    g = lambda c, w: np.exp(-0.5 * ((E - c) / w) ** 2)
    dTi = g(centres.get("Ti", -2.0), 0.6) * 4.0
    dV = g(centres.get("V", -1.0), 0.6) * 4.0
    pO = g(-5.0, 1.0)

    elements = []
    for el, c in zip(mx.species, mx.counts):
        elements += [el] * c

    os.makedirs(os.path.join(case, "DOS"), exist_ok=True)
    with open(os.path.join(case, "DOS", "DOSCAR"), "w") as f:
        f.write(f" {mx.nAtoms} {mx.nAtoms} 0 0\n 1000.0\n")
        for v in ("6 0 0", "-3 5.19 0", "0 0 27"):
            f.write(" " + v + "\n")
        f.write(" " + mx.label + "\n")
        f.write(f"  6.0 -10.0 {nedos} {EF}\n  energy dos int\n")
        tot = 8 * dTi + 4 * dV + 8 * pO
        for i in range(nedos):
            f.write(f"{Eabs[i]:11.5f}{tot[i]:13.5f}{tot[i] * 0.01:13.5f}\n")
        f.write("\n")
        for el in elements:
            if el == "Ti":
                d, s, p = dTi, 0.1 * dTi, 0.2 * dTi
            elif el == "V":
                d, s, p = dV, 0.1 * dV, 0.2 * dV
            elif el == "O":
                d, s, p = 0.01 * dV, 0.2 * pO, pO
            else:
                d, s, p = 0.02 * dV, 0.3 * pO, 0.5 * pO
            for i in range(nedos):
                f.write(f"{Eabs[i]:11.5f}{s[i]:13.5f}"
                        + "".join(f"{p[i] / 3:13.5f}" for _ in range(3))
                        + "".join(f"{d[i] / 5:13.5f}" for _ in range(5)) + "\n")
            f.write("\n")

    # The DOSCAR sits in <case>/DOS/ while the geometry lives under
    # <case root>/... -- copy it next to it so PDOS()'s upward search finds it.
    os.makedirs(case, exist_ok=True)
    shutil.copy(src, os.path.join(case, "POSCAR"))
    return os.path.join(case, "DOS", "DOSCAR"), src


def test_dband():
    print("\n[4] dband.py -- synthetic DOSCAR with a known eps_d")

    case = os.path.join(WORK, "synthetic_case")
    true_Ti, true_V = -2.0, -1.0
    dos, _ = make_synthetic(case, {"Ti": true_Ti, "V": true_V})
    d = PDOS(dos)

    check("atom list read from POSCAR, not from SYSTEM",
          d.elements.count("Ti") == 8 and d.elements.count("V") == 4)
    check("energies referenced to E_F", approx(d.E[-1], 6.0, 1e-6))

    got_Ti = d.dband_center("Ti")
    check(f"eps_d(Ti) = {true_Ti} on the occupied window",
          approx(got_Ti, true_Ti, 0.01), f"got {got_Ti:.4f}")
    check("eps_d(Ti) exact on the full-band window",
          approx(d.dband_center("Ti", emin=-10, emax=6), true_Ti, 5e-3))

    got_V = d.dband_center("V")
    check("eps_d(V) pulled toward E_F by the integration cutoff",
          true_V - 0.2 < got_V < true_V, f"got {got_V:.4f}")
    check("eps_d(V) lies above eps_d(Ti) (direction of the thesis claim)",
          got_V > got_Ti, f"V={got_V:.3f}  Ti={got_Ti:.3f}")

    m = d.d_moments("Ti")
    check("d-band width recovered (sigma = 0.6)", approx(m["width"], 0.6, 0.05),
          f"got {m['width']:.3f}")

    lay = d.layers(n=2)
    check("M sublattice split into 3 layers of 4",
          sorted(len(v) for v in lay.values()) == [4, 4, 4],
          str({k: len(v) for k, v in lay.items()}))

    s = d.summary()
    check("summary exposes layer-resolved eps_d",
          any(k.startswith("eps_d(") and "@" in k for k in s), str(list(s)))

    png = os.path.join(case, "pdos.png")
    d.plot(out=png)
    check("PDOS plot written", os.path.exists(png))


def test_parse_label():
    print("\n[5] dband.py -- label parsing")
    r = parse_label("Ti2.5V0.5C2O2_s2x2_sqs_seed42")
    check("host/dopant/y parsed",
          (r["host"], r["dopant"], r["y"]) == ("Ti", "V", 0.5), str(r))
    check("pattern + seed parsed", r["pattern"] == "sqs" and r["seed"] == 42, str(r))


# --------------------------------------------------------------------------- #
#  3. adsorption / neb / energetics
# --------------------------------------------------------------------------- #

def test_adsorption():
    print("\n[6] adsorption.py -- site selection and placement")
    from adsorption import Slab, make_adsorption_case, adsorption_energy, REFERENCES

    p_out = _slab_poscar("outer")
    p_in = _slab_poscar("inner")
    s_out, s_in = Slab(p_out), Slab(p_in)

    inv = s_out.site_inventory()
    check("8 surface O sites in a 2x2 cell", len(inv) == 8, str(len(inv)))
    check("outer pattern: every surface O sits above Ti",
          all(r["metal_below"] == "Ti" for r in inv),
          str({r["metal_below"] for r in inv}))
    check("inner pattern: every surface O sits above V",
          all(r["metal_below"] == "V" for r in s_in.site_inventory()),
          str({r["metal_below"] for r in s_in.site_inventory()}))

    ads = os.path.join(WORK, "ads.vasp")
    info = make_adsorption_case(p_out, ads, adsorbate="H", site="Otop",
                                both_sides=True)
    check("both_sides gives 2 adsorbates", info["n_adsorbates"] == 2, str(info))
    z = [p[2] for p in info["positions"]]
    check("mirrored adsorbate sits on the opposite face", max(z) - min(z) > 5.0,
          f"dz = {max(z) - min(z):.2f}")

    st = read_poscar(ads)
    check("adsorbate added as a new species", "H" in st["species"], str(st["species"]))
    check("atom count grew by 2", st["nAtoms"] == s_out.nAtoms + 2, str(st["nAtoms"]))

    REFERENCES["H"] = -3.4                     # placeholder: 1/2 E(H2)
    e = adsorption_energy(-500.0, -495.0, element="H", n_ads=2)
    check("E_ads arithmetic (per adsorbate)",
          approx(e, (-500.0 + 495.0 + 2 * 3.4) / 2, 1e-9), f"got {e}")
    try:
        adsorption_energy(-1, -2, element="Li")          # no reference set
        check("missing reference raises", False)
    except ValueError:
        check("missing reference raises a helpful error", True)


def test_neb():
    print("\n[7] neb.py -- image interpolation")
    from neb import make_neb
    from adsorption import make_adsorption_case

    slab = _slab_poscar("outer")
    ini = os.path.join(WORK, "neb_ini.vasp")
    fin = os.path.join(WORK, "neb_fin.vasp")
    make_adsorption_case(slab, ini, adsorbate="H", site="Otop",
                         o_index=24, both_sides=False)
    make_adsorption_case(slab, fin, adsorbate="H", site="Otop",
                         o_index=26, both_sides=False)

    out = make_neb(ini, fin, out=os.path.join(WORK, "NEB"), n_images=5,
                   potcar=os.path.join(os.path.dirname(os.path.dirname(slab)),
                                       "POTCAR"))
    dirs = [d for d in sorted(os.listdir(out)) if d.isdigit()]
    check("5 moving images + 2 endpoints = 7 dirs", len(dirs) == 7, str(dirs))
    ninc = parse_incar(os.path.join(out, "INCAR"))
    check("IMAGES = 5", ninc.get("IMAGES") == "5")
    check("LCLIMB on", ninc.get("LCLIMB") == ".TRUE.")
    # the adsorbate adds an H species -> 5 species in the LDAU arrays
    check("LDAUU written in POSCAR order (auto-detected)",
          ninc.get("LDAUU") == "2.50  2.50  0.00  0.00  0.00",
          str(ninc.get("LDAUU")))

    a = read_poscar(os.path.join(out, "00", "POSCAR"))
    b = read_poscar(os.path.join(out, "06", "POSCAR"))
    check("endpoint 00 == initial",
          np.allclose(a["positions"][-1], read_poscar(ini)["positions"][-1]))
    check("endpoint 06 == final",
          np.allclose(b["positions"][-1], read_poscar(fin)["positions"][-1]))
    mid = read_poscar(os.path.join(out, "03", "POSCAR"))
    check("middle image is halfway",
          approx(mid["positions"][-1][0],
                 0.5 * (a["positions"][-1][0] + b["positions"][-1][0]), 1e-6))


def test_energetics():
    print("\n[8] energetics.py -- mixing / site preference")
    from energetics import mixing_energy, site_preference, rank_configurations

    check("E_mix(x=0) == 0", approx(mixing_energy(-100.0, -100.0, -90.0, 0.0), 0.0))
    check("E_mix(x=1) == 0", approx(mixing_energy(-90.0, -100.0, -90.0, 1.0), 0.0))
    check("E_mix negative for a favoured solid solution",
          mixing_energy(-96.0, -100.0, -90.0, 0.5) < 0,
          f"{mixing_energy(-96.0, -100.0, -90.0, 0.5):.3f}")
    check("site_preference sign convention (positive = outer preferred)",
          approx(site_preference(-100.0, -99.5, 1), 0.5))

    rows = [{"y": 1.0, "pattern": "outer", "E": -100.2, "label": "a"},
            {"y": 1.0, "pattern": "inner", "E": -100.0, "label": "b"},
            {"y": 1.5, "pattern": "outer", "E": -99.0, "label": "c"}]
    r = rank_configurations(rows)
    check("ranking picks the lowest-energy pattern",
          r[0]["pattern"] == "outer" and r[0]["dE_min"] == 0.0, str(r[:2]))
    check("dE reported relative to the minimum", approx(r[1]["dE_min"], 0.2))


# --------------------------------------------------------------------------- #
#  4b. End-to-end: fabricated results tree -> full analysis pipeline
# --------------------------------------------------------------------------- #

def test_e2e_analyze():
    print("\n[4b] end-to-end analyse on a fabricated results tree")

    from dband import scan as dband_scan, plot_trend
    from energetics import collect_energies, rank_configurations, hull_table
    from make_fake_results import total_energy as fE, write_fake_outcar
    from make_fake_results import write_synthetic_doscar
    from alloy import AlloyMX as _A
    import numpy as np

    for y in (0.0, 0.5, 1.0, 1.5, 3.0):
        for pat in ("outer", "inner", "sqs"):
            if y == 1.5 and pat == "inner":
                continue                       # physically impossible
            mx = _A(y=y, pattern=pat, root=WORK)
            try:
                mx.write_inputs(encut=550, vdw="d3")
            except ValueError:
                continue
            case = os.path.join(WORK, mx.label, "ABC", "HM")
            write_fake_outcar(os.path.join(case, "opt", "OUTCAR"), fE(y, pat))
            with open(os.path.join(case, "opt", "POSCAR")) as f, \
                    open(os.path.join(case, "CONTCAR"), "w") as g:
                g.write(f.read())
            write_synthetic_doscar(case, y, pat, mx)

    rows = dband_scan(WORK)
    check("dband.scan finds every fabricated case (>= 12)",
          len(rows) >= 12, f"got {len(rows)}")

    erows = rank_configurations(collect_energies(WORK))
    check("rank_configurations runs on the fabricated tree",
          len(erows) >= 12, f"got {len(erows)}")
    erows = hull_table(erows, -600.0, -607.5, n=2)
    for y in (1.0, 1.5):
        grp = sorted([r for r in erows if r["y"] == y],
                     key=lambda r: r["E"])
        ps = [r["pattern"] for r in grp]
        check(f"y={y}: pattern ranking respects the model ordering",
              ps[0] == "outer" and "sqs" in ps[:2], str(ps))

    # plot_trend with extra series must not crash on numpy 2.x
    p = os.path.join(WORK, "trend.png")
    plot_trend(rows, key="eps_d(Ti)",
               extra={"C_g": [r.get("eps_d(Ti)", np.nan) for r in rows]},
               out=p)
    check("plot_trend runs with extra series on numpy 2.x",
          os.path.exists(p))


# --------------------------------------------------------------------------- #

def main(verbose: bool = False):
    os.makedirs(WORK, exist_ok=True)
    try:
        test_alloy_basics()
        test_alloy_symmetry()
        test_alloy_inputs()
        test_dband()
        test_parse_label()
        test_adsorption()
        test_neb()
        test_energetics()
        test_e2e_analyze()
    finally:
        if not verbose:
            shutil.rmtree(WORK, ignore_errors=True)
        else:
            print(f"\n(kept workdir: {os.path.join(ROOT, WORK)})")

    print("\n" + "=" * 62)
    print(f"  {len(PASSED)} passed, {len(FAILED)} failed")
    for f in FAILED:
        print(f"   FAILED: {f}")
    print("=" * 62)
    return 1 if FAILED else 0


if __name__ == "__main__":
    sys.exit(main(verbose="-v" in sys.argv))
