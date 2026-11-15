#!/usr/bin/env python3
"""Stand-in for the Vina CLI used by the test-suite (mimics argument handling/output)."""
import sys, time, os
a = sys.argv[1:]
if "--version" in a:
    print("AutoDock Vina v1.2.7"); sys.exit(0)
def opt(n):
    return a[a.index(n) + 1] if n in a else None
cfg = {}
for line in open(opt("--config")):
    if "=" in line and not line.startswith("#"):
        k, v = line.split("=", 1); cfg[k.strip()] = v.strip()
for req in ("center_x", "size_z", "exhaustiveness", "num_modes", "scoring"):
    if req not in cfg:
        print("Error: missing %s" % req); sys.exit(1)
lig = open(opt("--ligand")).read()
if "BADLIG" in lig:
    print("Parse error on line 3 in file \"%s\": Unknown AutoDock atom type" % opt("--ligand")); sys.exit(1)
import random
seed = int(opt("--seed")) if opt("--seed") else random.randint(-2**31, 2**31 - 1)
print("Performing docking (random seed: %d) ... " % seed)
sys.stdout.write("0%   10   20   30   40   50   60   70   80   90   100%\n|----|----|----|----|----|----|----|----|----|----|\n")
for _ in range(51):
    sys.stdout.write("*"); sys.stdout.flush(); time.sleep(0.002)
sys.stdout.write("\n")
n = int(cfg["num_modes"]); base = -5.0 - (len(lig) % 40) / 10.0
with open(opt("--out"), "w") as f:
    for i in range(min(n, 3)):
        f.write("MODEL %d\nREMARK VINA RESULT:   %6.3f      0.000      0.000\n" % (i + 1, base + 0.4 * i))
        f.write("ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00    +0.000 C \nENDMDL\n")
print("mode |   affinity | dist from best mode\n     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n-----+------------+----------+----------")
