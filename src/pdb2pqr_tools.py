"""PDB2PQR + PROPKA: pKa-based protonation states, run in the private preparation environment.

Pure Python (no ChimeraX imports). PDB2PQR is called as ``python -m pdb2pqr``; the pKa summary is read from its
output, and the residue names it writes (HID/HIE/HIP, ASH, GLH, LYN, CYM ...) are compared with the input.
"""

import csv
import os
import re

from .prep_env import PrepError, run_stream

CLOSE = 1.0     # |pKa - pH| <= 1: strong flag (PROPKA is accurate to about one pKa unit)
NEAR = 2.0      # |pKa - pH| <= 2: soft flag, for extra care
ACIDIC = ("ASP", "GLU", "CYS", "TYR", "C-")

_ROW = re.compile(r"^\s*([A-Z][A-Z0-9+\-]{1,3})\s+(-?\d+)\s+(\S)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)")


def pdb2pqr_args(ph, pdb_out, pqr_out, in_pdb, ff="AMBER"):
    """Arguments after ``python -m pdb2pqr``. AMBER naming gives HID/HIE/HIP, ASH, GLH, LYN, CYM."""
    return ["--ff", ff, "--ffout", ff, "--titration-state-method", "propka", "--with-ph", "%g" % ph,
            "--keep-chain", "--pdb-output", pdb_out, in_pdb, pqr_out]


def classify(pka, ph):
    """'close' (within 1 pH unit), 'near' (within 2) or ''."""
    d = round(abs(pka - ph), 6)      # rounding: 7.4 - 6.4 must count as exactly 1
    return "close" if d <= CLOSE else ("near" if d <= NEAR else "")


def parse_pka_summary(text):
    """Rows of the 'SUMMARY OF THIS PREDICTION' table: dicts with group, number, chain, pka, model_pka."""
    rows, inside = [], False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            inside, rows = True, []
            continue
        if not inside:
            continue
        m = _ROW.match(line)
        if m:
            rows.append({"group": m.group(1), "number": int(m.group(2)), "chain": m.group(3),
                         "pka": float(m.group(4)), "model_pka": float(m.group(5))})
        elif rows and not line.strip():
            break
    return rows


def residue_names(pdb_path):
    """{(chain, number): residue name} of the ATOM/HETATM records of a PDB file."""
    out = {}
    try:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
                    out.setdefault((line[21].strip() or " ", int(line[22:26])), line[17:20].strip())
    except (OSError, ValueError):
        pass
    return out


def strip_all_hydrogens(pdb_in, pdb_out):
    """Copy a PDB file without hydrogen atoms (PDB2PQR adds its own). Returns the number removed."""
    n = 0
    with open(pdb_in, encoding="utf-8", errors="replace") as f, open(pdb_out, "w", encoding="utf-8") as g:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
                elem = line[76:78].strip() if len(line) >= 78 else ""
                name = line[12:16].strip().lstrip("0123456789")
                if elem in ("H", "D") or (not elem and line.startswith("ATOM") and name[:1] == "H"):
                    n += 1
                    continue
            g.write(line)
    return n


def copy_seqres(src_pdb, dst_pdb):
    """Prepend the SEQRES records of *src_pdb* to *dst_pdb* (PDB2PQR drops all header records). ChimeraX uses
    SEQRES to tell real chain ends from cut ones. Returns the number of records copied."""
    try:
        with open(src_pdb, encoding="utf-8", errors="replace") as f:
            seq = [l for l in f if l.startswith("SEQRES")]
        if not seq:
            return 0
        with open(dst_pdb, encoding="utf-8", errors="replace") as f:
            body = f.read()
        if "\nSEQRES" in "\n" + body:
            return 0
        with open(dst_pdb, "w", encoding="utf-8") as f:
            f.write("".join(seq) + body)
        return len(seq)
    except OSError:
        return 0


def build_table(rows, before, after, ph):
    """Add flag, predicted charge state and the name PDB2PQR applied to each pKa row."""
    table = []
    for r in rows:
        key = (r["chain"], r["number"])
        old, new = before.get(key), after.get(key)
        acidic = r["group"] in ACIDIC
        deprot = r["pka"] < ph
        predicted = "deprotonated" if deprot else "protonated"
        if r["group"] in ("N+", "C-"):
            label = "N-terminus" if r["group"] == "N+" else "C-terminus"
            residue = "%s %s%d" % (label, r["chain"], r["number"])
        else:
            residue = "%s %s%d" % (r["group"], r["chain"], r["number"])
        row = dict(r)
        row.update({"residue": residue, "delta": r["pka"] - ph, "flag": classify(r["pka"], ph),
                    "predicted": predicted, "acidic": acidic,
                    "renamed": ("%s → %s" % (old, new)) if old and new and old != new else ""})
        table.append(row)
    return table


def write_pka_csv(path, table, ph):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["residue", "pKa", "model_pKa", "pKa_minus_pH", "pH", "predicted", "applied_name_change", "check"])
        for r in table:
            w.writerow([r["residue"], "%.2f" % r["pka"], "%.2f" % r["model_pka"], "%+.2f" % r["delta"], "%g" % ph,
                        r["predicted"], r["renamed"],
                        {"close": "within 1 pH unit", "near": "within 2 pH units"}.get(r["flag"], "")])


def run_pdb2pqr(env, in_pdb, out_dir, base, ph, log=None, cancelled=None, timeout=1800):
    """Run PDB2PQR/PROPKA on *in_pdb* (hydrogens are removed first). Returns a dict:
    pdb (protonated PDB), table (see build_table), renamed (list), log_path, csv_path, log_text.
    """
    os.makedirs(out_dir, exist_ok=True)
    heavy = os.path.join(out_dir, base + "_input.pdb")
    out_pdb = os.path.join(out_dir, base + "_pdb2pqr.pdb")
    out_pqr = os.path.join(out_dir, base + "_pdb2pqr.pqr")
    for p in (out_pdb, out_pqr):
        if os.path.exists(p):
            os.remove(p)
    removed = strip_all_hydrogens(in_pdb, heavy)
    if log:
        log("Running PDB2PQR + PROPKA at pH %g (%d existing hydrogens removed first) ..." % (ph, removed))
    lines = []

    def on_line(line):
        lines.append(line)

    cmd = env.python_command(["-m", "pdb2pqr"] + pdb2pqr_args(ph, out_pdb, out_pqr, heavy))
    code = run_stream(cmd, on_line, cancelled, cwd=out_dir, env=env.environ(), timeout=timeout)
    text = "\n".join(lines)
    log_path = os.path.join(out_dir, base + "_pdb2pqr.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("pdb2pqr %s\n\n%s\n" % (" ".join(cmd[-12:]), text))
    if code == -1:
        raise PrepError("Cancelled.")
    if code != 0 or not os.path.isfile(out_pdb):
        errs = [l for l in lines if "ERROR" in l or "CRITICAL" in l]
        raise PrepError("PDB2PQR failed (exit code %d)%s. Full log: %s" % (
            code, ": " + errs[-1].split(":", 3)[-1].strip()[:200] if errs else "", log_path))
    rows = parse_pka_summary(text)
    table = build_table(rows, residue_names(heavy), residue_names(out_pdb), ph)
    csv_path = os.path.join(out_dir, base + "_pka.csv")
    write_pka_csv(csv_path, table, ph)
    renamed = [r["residue"] + ": " + r["renamed"] for r in table if r["renamed"]]
    return {"pdb": out_pdb, "table": table, "renamed": renamed, "log_path": log_path, "csv_path": csv_path,
            "ph": ph, "log_text": text, "removed_h": removed}
