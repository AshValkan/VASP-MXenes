"""
run_alloy.py -- config-driven entry point for the alloy workflow.

Reads ``alloy_config.yaml`` and does everything that does not need a
supercomputer: builds the case matrix, writes the VASP inputs, and (later,
once the runs are finished) drives the analysis.

    python run_alloy.py --config alloy_config.yaml            # generate inputs
    python run_alloy.py --config alloy_config.yaml --analyze  # analyse results
    python run_alloy.py --config alloy_config.yaml --dry-run  # list, don't write

Why a config file
-----------------
Every number in the thesis has to be traceable.  With everything in one
YAML, "which U did you use for V?" has a one-line answer, and re-running the
whole study after a change is a one-liner instead of an archaeology exercise
through half a dozen scripts.
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def load_config(path: str) -> dict:
    try:
        import yaml
    except ImportError:
        sys.exit("PyYAML is required for the config file:  pip install pyyaml")
    with open(path) as f:
        return yaml.safe_load(f)


# --------------------------------------------------------------------------- #

def cmd_generate(cfg: dict, dry: bool = False) -> None:
    from alloy import AlloyMX, generate_batch, suggest_supercell

    s = cfg["system"]
    ys = list(s["ys"])
    sc = s.get("supercell", [2, 2, 1])

    rec = suggest_supercell(ys)
    print("=" * 66)
    print(f"system      : {s['host']}_(3-y){s['dopant']}_y {s['X']}{s['n']} "
          f"{s['T'] or '(pristine)'}   {s['stacking']}-{s['hollows']}")
    print(f"compositions: y = {ys}")
    print(f"patterns    : {s['patterns']}")
    print(f"supercell   : {sc[0]}x{sc[1]}x{sc[2]}  "
          f"({(s['n'] + 1) * sc[0] * sc[1]} M sites)")
    print(f"smallest valid cell for these y: {rec}   "
          f"(you asked for {sc[0]}x{sc[1]})")
    if rec and rec != tuple(sc[:2]):
        print("  note: a smaller cell is commensurate; using the one you "
              "specified.")
    d = cfg["dft"]
    print(f"DFT         : ENCUT={d['encut']}  vdW={d['vdw']}  "
          f"U={d.get('hubbard_u', {})}  ISPIN={2 if d.get('spin', True) else 1}")
    print("=" * 66)

    if dry:
        print("\n[dry run] nothing written.  "
              f"Would create {len(ys) * len(s['patterns'])} cases under "
              f"{cfg['output']['root']}/")
        return

    made = generate_batch(
        ys, s["patterns"], host=s["host"], dopant=s["dopant"], n=s["n"],
        X=s["X"], T=s["T"], supercell=tuple(sc), seed=s.get("seed", 42),
        stacking=s["stacking"], hollows=s["hollows"],
        v=cfg["geometry"].get("vacuum", 20.0), root=cfg["output"]["root"],
        encut=d["encut"], vdw=d["vdw"],
        spin=d.get("spin", True), ediff=d.get("ediff", "1E-06"),
        ediffg=d.get("ediffg", "-0.01"), sigma=d.get("sigma", 0.05),
        lorbit=d.get("lorbit", 11), nedos=d.get("nedos", 5000),
        dipol=d.get("dipole", True))

    # Per-case U overrides (e.g. a convergence study in U).
    us = d.get("hubbard_u", {})
    if us:
        for mx in made:
            mx.set_U(**us)

    print(f"\n{len(made)} cases written under {cfg['output']['root']}/")
    print("Next: copy the tree to the cluster, put your own queue script in "
          "car/script,")
    print("      then run opt.py / calculate.py as upstream describes.")


def cmd_analyze(cfg: dict, dry: bool = False) -> None:
    from analyze_alloy import run
    o = cfg["output"]
    d = cfg.get("dband", {})
    run(root=o["root"], n=cfg["system"]["n"],
        m_elements=d.get("elements", ("Ti", "V")),
        emin=d.get("emin", -10.0), emax=d.get("emax", 0.0),
        exp=cfg.get("experimental", {}).get("capacitance_file"),
        ads=cfg.get("adsorption", {}).get("results_file"),
        outdir=o.get("results", "results"), dpi=o.get("dpi", 300))


# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="alloy_config.yaml")
    ap.add_argument("--analyze", action="store_true",
                    help="analyse a finished results tree instead of generating inputs")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if not os.path.exists(a.config):
        sys.exit(f"config not found: {a.config}")
    cfg = load_config(a.config)

    if a.analyze:
        cmd_analyze(cfg, dry=a.dry_run)
    else:
        cmd_generate(cfg, dry=a.dry_run)


if __name__ == "__main__":
    main()
