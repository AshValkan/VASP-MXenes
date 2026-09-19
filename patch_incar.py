import os, re, glob

base = os.path.expanduser("~/calc/VASP-MXenes-alloy/MXenes_alloy")
cases = [
    "Ti3V0C2O2_s2x2_outer",
    "Ti2.5V0.5C2O2_s2x2_outer",
    "Ti2V1C2O2_s2x2_outer",
    "Ti1.5V1.5C2O2_s2x2_outer",
]

# (regex to match a line, replacement), applied in order
subs = [
    (re.compile(r"^EDIFF\s*=.*$", re.M), "EDIFF  = 1e-05"),
    (re.compile(r"^EDIFFG\s*=.*$", re.M), "EDIFFG = -0.02"),
    (re.compile(r"^NWRITE\s*=.*$", re.M), "NWRITE = 1"),
    (re.compile(r"^LDAUPRINT\s*=.*$", re.M), "LDAUPRINT = 1"),
]

for c in cases:
    inc = os.path.join(base, c, "ABC/HM/opt/INCAR")
    if not os.path.exists(inc):
        print("SKIP (no INCAR):", inc)
        continue
    with open(inc) as f:
        txt = f.read()
    for rx, rep in subs:
        txt, n = rx.subn(rep, txt)
        if n == 0:
            print("  WARN did not replace for", rep, "in", c)
    # insert AMIX/BMIX after NELM line if not already present
    if "AMIX" not in txt:
        txt = re.sub(r"(^NELM\s*=.*$\n)", r"\1AMIX   = 0.2\nBMIX   = 0.0001\n", txt, flags=re.M)
    with open(inc, "w") as f:
        f.write(txt)
    print("PATCHED:", inc)

print("=== verify ===")
for c in cases:
    inc = os.path.join(base, c, "ABC/HM/opt/INCAR")
    with open(inc) as f:
        for line in f:
            if re.match(r"^\s*(EDIFF|EDIFFG|NWRITE|LDAUPRINT|AMIX|BMIX|NELM)\s*=", line):
                print(c, ":", line.rstrip())
