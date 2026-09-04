# Alloyed MXene extension (fork of diegonti/VASP-MXenes)

> **Upstream README** (`README.md`) documents the original single-element
> workflow and is unchanged. This file documents what **this fork** adds.

This fork re-targets the upstream toolkit at the system studied in

> **合金化 MXene 电子结构调控与超级电容器储能性能研究**
> (Alloyed-MXene electronic-structure modulation for supercapacitor energy storage)
> Wang Yaoyu, NUIST, supervisor Prof. Yu Feng

i.e. **Ti<sub>3−y</sub>V<sub>y</sub>C<sub>2</sub>O<sub>2</sub>** — a *binary* M-sublattice
solid solution, in supercells, with the d-band centre as the central descriptor.

Upstream handles single-metal Mn+1XnT2 in the primitive cell, and stops at DOS /
band structure / work function. Everything needed for a composition series,
site-resolved electronic structure and adsorption/diffusion is added here.

<br>

---

## 1. What was added

| file | purpose | upstream equivalent |
|---|---|---|
| `alloy.py` | binary-M MXenes, supercells, dopant site control | `VASP.py` (extended) |
| `dband.py` | **d-band centre**, band moments, composition trends | `DOS.py` (replaced for this use) |
| `adsorption.py` | H/Li/Na/K placement + E_ads | — (new) |
| `neb.py` | CI-NEB setup + barrier extraction | — (new) |
| `energetics.py` | E_form, E_mix, site-preference energy | — (new) |
| `analyze_alloy.py` | one command → thesis figures | `analyzer.py` (superseded) |
| `run_alloy.py` | config-driven entry point | — (new) |
| `alloy_config.yaml` | the single source of truth for the system | — (new) |
| `tests/selftest.py` | 63 checks, runs **without VASP** | — (new) |

<br>

## 2. Quick start

```bash
pip install numpy matplotlib pyyaml        # pymatgen optional (IDPP for NEB)

python run_alloy.py --config alloy_config.yaml             # write all inputs
python run_alloy.py --config alloy_config.yaml --dry-run   # just list them
python tests/selftest.py                                   # 63 checks, no VASP
```

This creates, for every composition × pattern:

```
MXenes_alloy/
└── Ti2V1C2O2_s2x2_outer/
    └── ABC/HM/
        ├── opt/    POSCAR POTCAR KPOINTS INCAR     # relax
        ├── DOS/    ...  LORBIT=11, NEDOS=5000      # -> dband.py
        ├── BS/     ...
        ├── WF/     LVTOT, LVHAR, dipole corr.      # -> work function
        └── bader/  LAECHG                          # -> Bader charges
```

Then the upstream cluster workflow applies unchanged (`opt.py`,
`calculate.py`), and the new analysis side is:

```bash
python dband.py --scan MXenes_alloy --plot        # eps_d vs composition
python analyze_alloy.py --root MXenes_alloy \
       --exp measured.json --ads adsorption.json  # thesis figures
python energetics.py --root MXenes_alloy \
       --host-end -123.45 --dopant-end -118.90    # E_mix / ranking
```

<br>

## 3. The three questions this fork is built to answer

**Q1 — where does the V actually go?**

For n = 2 the M sublattice is *outer / inner / outer*, and the two are not
equivalent. `pattern=outer|inner|sqs|...` generates them explicitly, and
`energetics.site_preference()` turns the comparison into a number:

```python
from energetics import site_preference
dE = site_preference(E_outer_case, E_inner_case, n_dopant=4)
#   dE > 0 -> V prefers the surface layer
#   dE < 0 -> V prefers the inner layer
```

> ⚠️ **Physics note that will save you a week.** In an H_M (FCC) termination
> the surface oxygen sits **directly above the middle metal layer**. So
> `pattern=inner` puts V *under every adsorption site*, while `pattern=outer`
> leaves the adsorption sites on Ti. The two patterns change adsorption
> through completely different mechanisms — do not treat them as two random
> guesses, and `adsorption.Slab.site_inventory()` prints exactly this
> metal-under-each-site table.

**Q2 — does the d-band centre actually move toward E_F?**

`dband.PDOS.dband_center()` gives ε_d per element *and per layer*, and
`plot_trend()` produces the composition plot with a linear fit and **R²**.
Without that R², "V moves ε_d up ⇒ capacitance improves" is a picture, not a
result.

The integration window defaults to the occupied band (−10 eV … E_F), the
standard convention. Whatever you pick, keep it identical across every
composition and state it in the thesis.

**Q3 — how does that become rate capability?**

`neb.py` builds CI-NEB runs for the proton/cation hop (needs the VTST binary),
and `barrier_from_neb()` extracts E_a forward and backward. Capacity and rate
pull in opposite directions; showing V improves both is the interesting part.

<br>

## 4. Bugs fixed relative to upstream

These only show up once M is not a single element, but they are silent when
they do — worth knowing if you ever merge upstream changes.

1. **POTCAR selection** — `VASP.MX.POTCAR` special-cases only
   `self.atoms[0]`, so a second metal silently gets the *bare* PAW potential
   (`V` instead of `V_pv`). Fixed: every species is mapped through
   `alloy.PP_ALIAS`.
2. **`n` inference** — `splitNums("Ti2V1C2O2")` makes `index[0] == 2`, so
   `self.n` would be 1, not 2. `AlloyMX` takes `n` explicitly.
3. **No vdW / no Hubbard U** — layered MXenes need DFT-D3 and the thesis
   needs U<sub>eff</sub> = 2.5 eV on both Ti-3d and V-3d; the LDAU arrays are
   now written **in POSCAR species order** (getting that wrong is silent in
   VASP).
4. **KPOINTS** — a supercell kept the primitive mesh, i.e. ~4× the cost.
   Now divided by sx, sy.
5. **Direct/Cartesian mislabelled** — the new adsorption/NEB writers
   originally emitted fractional coordinates under a `Cartesian` header, which
   VASP reads as Cartesian without complaining. Caught by
   `tests/selftest.py`; both now write `Direct`.

<br>

## 5. Things that are deliberately *not* automatic

* **Reference energies** (`adsorption.REFERENCES`) are placeholders. Compute
  ½E(H₂) and the bulk metals *with your own pseudopotentials*; quoting a
  reference from someone else's paper is not defensible in a thesis.
* **The SQS here is a heuristic** (best-of-N random draws scored by
  pair-disorder). For 12–36 sites that is entirely sufficient, but if a
  reviewer asks for a rigorous SQS, generate it with ATAT `mcsqs` or ICET and
  drop the POSCAR in place — the rest of the pipeline does not care where the
  POSCAR came from.
* **V₃C₂Tₓ cannot be synthesised.** It is used here only as the computational
  end member for E_mix. Say so in the thesis — an MXene reviewer will know.
* **Queue scripts** (`car/script*`) are upstream's. Put your own cluster's
  there.

<br>

## 6. Running the tests

```bash
python tests/selftest.py      # 63 checks, ~5 s, no VASP needed
```

Covers: composition/atom-count invariants, surface symmetry (no spurious slab
dipole), POTCAR species+order, INCAR contents (ENCUT/IVDW/LDAUU/LDAUL/MAGMOM/
LORBIT/dipole), k-mesh reduction, DOSCAR parsing against an analytically known
ε_d, layer resolution, adsorption site selection by underlying metal, NEB
image construction and endpoint preservation, and the mixing-energy algebra.

<br>

## 7. Repo layout notes

```
upstream  = https://github.com/diegonti/VASP-MXenes   (fetch/push)
branch    = alloy-mxene-tiv
```

To publish your own fork:

```bash
# on GitHub: fork diegonti/VASP-MXenes, then
git remote add origin https://github.com/<YOU>/VASP-MXenes.git
git push -u origin alloy-mxene-tiv
```

<br>

---

*Fork extensions: Wang Yaoyu (王耀宇), Nanjing University of Information Science
and Technology, 2026. Upstream toolkit: Diego Ontiveros —
[github.com/diegonti/VASP-MXenes](https://github.com/diegonti/VASP-MXenes).*
