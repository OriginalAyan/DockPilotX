"""PDBQT inspection, ligand discovery and multi-ligand library splitting.

Pure Python (no ChimeraX imports).
"""

import os
import re
import shutil

import numpy as np

# AutoDock atom types accepted by Vina (case as written by Meeko/ADFR/ADT).
VALID_AD_TYPES = {
    "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F", "MG", "MN", "P", "S", "SA",
    "CL", "CA", "FE", "ZN", "BR", "I", "Z", "G", "GA", "J", "Q", "CG0", "CG1", "CG2", "CG3",
    "G0", "G1", "G2", "G3", "W", "B", "SI", "AS", "AU", "AL", "CU", "CO", "NI", "K", "NA+",
}

_NAME_RE = re.compile(r"^REMARK\s+Name\s*=\s*(.+?)\s*$", re.IGNORECASE)
_SAFE_RE = re.compile(r"[^A-Za-z0-9._+-]+")


def _rec(line):
    """Record name (first whitespace-delimited token) of a PDBQT line."""
    parts = line.split(None, 1)
    return parts[0] if parts else ""


def safe_name(text, fallback="ligand"):
    """Make a string safe to use as a file stem on every platform."""
    s = _SAFE_RE.sub("_", (text or "").strip()).strip("._")
    return s[:80] if s else fallback


def _atom_type(line):
    """AutoDock atom type from an ATOM/HETATM record (columns 78-79)."""
    t = line[77:79].strip() if len(line) >= 78 else ""
    if not t:
        parts = line.split()
        t = parts[-1] if parts else ""
    return t


def inspect_pdbqt(path, max_bytes=None):
    """Summarize a PDBQT file without loading it all (streams line by line).

    Returns a dict: n_atoms, n_hydrogens, atom_types (set), n_root, n_torsdof,
    n_models, has_branch.
    """
    info = {"n_atoms": 0, "n_hydrogens": 0, "atom_types": set(), "n_root": 0,
            "n_torsdof": 0, "n_models": 0, "has_branch": False}
    with open(path, "r", errors="replace") as f:
        for line in f:
            rec = _rec(line)
            if rec in ("ATOM", "HETATM"):
                info["n_atoms"] += 1
                t = _atom_type(line)
                info["atom_types"].add(t)
                if t.upper() in ("H", "HD", "HS"):
                    info["n_hydrogens"] += 1
            elif rec == "ROOT":
                info["n_root"] += 1
            elif rec == "TORSDOF":
                info["n_torsdof"] += 1
            elif rec == "MODEL":
                info["n_models"] += 1
            elif rec == "BRANCH":
                info["has_branch"] = True
    return info


def check_receptor_pdbqt(path):
    """Return a list of (level, message); level is 'error' or 'warning'."""
    msgs = []
    if not os.path.isfile(path):
        return [("error", "Receptor file not found: %s" % path)]
    if not path.lower().endswith(".pdbqt"):
        msgs.append(("warning", "Receptor file does not end in .pdbqt."))
    info = inspect_pdbqt(path)
    if info["n_atoms"] == 0:
        return [("error", "Receptor has no ATOM/HETATM records - not a valid PDBQT.")]
    unknown = sorted(t for t in info["atom_types"] if t.upper() not in VALID_AD_TYPES)
    if unknown:
        msgs.append(("error", "Receptor has atom types Vina does not recognise: %s. "
                              "It was probably not written by a PDBQT preparation tool "
                              "(Meeko, ADFR, AutoDockTools)." % ", ".join(unknown)))
    if info["n_root"] or info["n_torsdof"]:
        msgs.append(("warning", "Receptor contains ROOT/TORSDOF records; it looks like a flexible "
                                "ligand or a receptor with flexible residues. A normal rigid "
                                "receptor PDBQT has neither."))
    if info["n_models"] > 1:
        msgs.append(("error", "Receptor has several MODEL records; supply a single rigid structure."))
    if info["n_hydrogens"] == 0:
        msgs.append(("warning", "Receptor has no polar hydrogens (HD). Hydrogen-bond scoring will be "
                                "poor; prepare the receptor with hydrogens added."))
    return msgs


def check_ligand_block(lines):
    """Validate one ligand's PDBQT lines. Returns a list of error strings."""
    errs = []
    atoms = [l for l in lines if _rec(l) in ("ATOM", "HETATM")]
    if not atoms:
        errs.append("no ATOM/HETATM records")
    if not any(l.startswith("ROOT") for l in lines):
        errs.append("no ROOT record (not a ligand PDBQT)")
    if not any(l.startswith("TORSDOF") for l in lines):
        errs.append("no TORSDOF record (file truncated?)")
    unknown = sorted({_atom_type(l) for l in atoms if _atom_type(l).upper() not in VALID_AD_TYPES})
    if unknown:
        errs.append("unrecognised atom types: %s" % ", ".join(unknown))
    return errs


def iter_pdbqt_blocks(fileobj):
    """Yield (name_or_None, [lines]) for every ligand in a (multi-)ligand PDBQT.

    A ligand ends at its TORSDOF record. MODEL/ENDMDL wrapper lines are
    dropped so that each yielded block is a self-contained ligand file.
    """
    block, name = [], None
    for raw in fileobj:
        line = raw.rstrip("\r\n")
        rec = _rec(line)
        if rec in ("MODEL", "ENDMDL"):
            continue
        if not line.strip() and not block:
            continue
        block.append(line)
        if name is None:
            m = _NAME_RE.match(line)
            if m:
                name = m.group(1)
        if rec == "TORSDOF":
            yield name, block
            block, name = [], None
    if any(_rec(l) in ("ATOM", "HETATM") for l in block):
        yield name, block  # trailing block without TORSDOF; reported later by check_ligand_block


def count_ligands(path):
    """Count ligands in a PDBQT file = number of TORSDOF records (streams the file)."""
    n = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b"TORSDOF"):
                n += 1
    return n


HBA_TYPES = {"OA", "NA", "SA"}


def _xyz(line):
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except ValueError:
        parts = line.split()
        nums = []
        for p in parts[5:]:
            try:
                nums.append(float(p))
            except ValueError:
                continue
        return tuple(nums[1:4]) if len(nums) >= 4 else (0.0, 0.0, 0.0)


def _charge(line):
    try:
        return float(line[70:76])
    except ValueError:
        parts = line.split()
        try:
            return float(parts[-2])
        except (ValueError, IndexError):
            return 0.0


def ligand_properties(lines):
    """Properties derivable from one ligand's PDBQT lines.

    Dimensions are the extents of the atom centres along the molecule's own
    principal axes (so they do not depend on how the file is oriented),
    sorted largest first. 'longest' is the largest atom-to-atom distance.
    These describe the conformation stored in the file; a flexible ligand
    can be longer when its torsions are extended.
    """
    atoms = [l for l in lines if _rec(l) in ("ATOM", "HETATM")]
    xyz = np.array([_xyz(l) for l in atoms], dtype=float) if atoms else np.zeros((0, 3))
    types = [_atom_type(l).upper() for l in atoms]
    heavy = [i for i, t in enumerate(types) if t not in ("H", "HD", "HS")]
    props = {
        "n_atoms": len(heavy),
        "n_polar_h": sum(1 for t in types if t == "HD"),
        "n_acceptors": sum(1 for t in types if t in HBA_TYPES),
        "charge": round(sum(_charge(l) for l in atoms), 3),
        "torsions": None, "n_branches": sum(1 for l in lines if _rec(l) == "BRANCH"),
        "dims": (0.0, 0.0, 0.0), "aabb": (0.0, 0.0, 0.0), "longest": 0.0, "smiles": "",
    }
    for l in lines:
        r = _rec(l)
        if r == "TORSDOF":
            try:
                props["torsions"] = int(l.split()[1])
            except (IndexError, ValueError):
                pass
        elif r == "REMARK" and "SMILES" in l and "IDX" not in l and not props["smiles"]:
            parts = l.split()
            if len(parts) >= 3 and parts[1] == "SMILES":
                props["smiles"] = parts[2]
    if len(heavy) >= 1:
        pts = xyz[heavy]
        props["aabb"] = tuple(float(v) for v in (pts.max(axis=0) - pts.min(axis=0)))
        if len(pts) >= 2:
            c = pts - pts.mean(axis=0)
            _u, _s, vt = np.linalg.svd(c, full_matrices=False)
            proj = c @ vt.T
            ext = proj.max(axis=0) - proj.min(axis=0)
            ext = np.sort(np.pad(ext, (0, 3 - len(ext))))[::-1]
            props["dims"] = tuple(float(v) for v in ext)
            diff = pts[:, None, :] - pts[None, :, :]
            props["longest"] = float(np.sqrt((diff ** 2).sum(axis=2)).max())
    return props


class LigandRecord:
    """One ligand found in an input file (a single-ligand file or one entry of a library)."""

    def __init__(self, name, source, index, single, props, problems):
        self.name = name
        self.source = os.path.abspath(source)
        self.index = index          # position of the ligand inside its file (0-based)
        self.single = single        # True: the whole file is this ligand
        self.props = props
        self.problems = problems

    @property
    def key(self):
        return (os.path.normcase(self.source), self.index)

    def __repr__(self):
        return "LigandRecord(%r, %r, #%d)" % (self.name, self.source, self.index)


class LigandItem:
    """A staged ligand file ready for docking."""

    def __init__(self, name, path, source=None, problems=None):
        self.name = name
        self.path = path
        self.source = source or path
        self.problems = problems or []

    def __repr__(self):
        return "LigandItem(%r, %r)" % (self.name, self.path)


def _unique(name, used):
    base, k, cand = name, 2, name
    while cand.lower() in used:
        cand = "%s_%d" % (base, k)
        k += 1
    used.add(cand.lower())
    return cand


def list_pdbqt_files(inputs, recursive=False):
    """Expand files/folders into a sorted list of .pdbqt file paths."""
    files = []
    for p in inputs:
        p = os.path.expanduser(str(p).strip().strip('"'))
        if os.path.isdir(p):
            if recursive:
                for root, _dirs, names in os.walk(p):
                    files += [os.path.join(root, n) for n in names if n.lower().endswith(".pdbqt")]
            else:
                files += [os.path.join(p, n) for n in os.listdir(p) if n.lower().endswith(".pdbqt")]
        elif os.path.isfile(p):
            files.append(p)
        else:
            raise FileNotFoundError("Ligand input not found: %s" % p)
    seen, out = set(), []
    for f in sorted(files, key=lambda s: s.lower()):
        key = os.path.normcase(os.path.abspath(f))
        if key not in seen:
            seen.add(key)
            out.append(os.path.abspath(f))
    return out


def _has_model_wrappers(path):
    with open(path, "r", errors="replace") as f:
        for line in f:
            if line[:5] in ("MODEL", "ENDMD"):
                return True
    return False


def scan_ligands(inputs, recursive=False, used=None, progress=None, skip_keys=()):
    """Read ligand PDBQT inputs and return a list of LigandRecord (no files are written).

    *used* is a set of lower-case names already taken (updated in place) so that
    names stay unique across several calls. Records whose key is in
    *skip_keys* (already listed by the caller) are not returned again.
    """
    used = used if used is not None else set()
    skip = set(skip_keys)
    files = list_pdbqt_files(inputs, recursive)
    records = []
    for i, src in enumerate(files):
        stem = safe_name(os.path.splitext(os.path.basename(src))[0], "ligand")
        n_lig = count_ligands(src)
        with open(src, "r", errors="replace") as f:
            blocks = list(iter_pdbqt_blocks(f))
        if not blocks:
            rec = LigandRecord(_unique(stem, used), src, 0, True, ligand_properties([]),
                               ["no ligand records (ATOM/ROOT/TORSDOF) found"])
            if rec.key not in skip:
                records.append(rec)
        elif len(blocks) == 1 and n_lig <= 1 and not _has_model_wrappers(src):
            lines = blocks[0][1]
            rec = LigandRecord(None, src, 0, True, ligand_properties(lines), check_ligand_block(lines))
            if rec.key not in skip:
                rec.name = _unique(stem, used)
                records.append(rec)
        else:
            for idx, (lname, lines) in enumerate(blocks):
                base = safe_name(lname, "%s_%d" % (stem, idx + 1)) if lname else "%s_%d" % (stem, idx + 1)
                rec = LigandRecord(None, src, idx, False, ligand_properties(lines), check_ligand_block(lines))
                if rec.key in skip:
                    continue
                rec.name = _unique(base, used)
                records.append(rec)
        if progress:
            progress(i + 1, len(files))
    return records


def stage_records(records, dest_dir):
    """Write one PDBQT per record into dest_dir (single files are copied, libraries streamed).

    Returns a list of LigandItem in the order of *records*.
    """
    os.makedirs(dest_dir, exist_ok=True)
    by_source = {}
    for r in records:
        by_source.setdefault(r.source, []).append(r)
    items = {}
    for src, recs in by_source.items():
        wanted = {r.index: r for r in recs}
        if all(r.single for r in recs):
            r = recs[0]
            dst = os.path.join(dest_dir, r.name + ".pdbqt")
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copyfile(src, dst)
            items[id(r)] = LigandItem(r.name, dst, src, r.problems)
            continue
        with open(src, "r", errors="replace") as f:
            for idx, (_lname, lines) in enumerate(iter_pdbqt_blocks(f)):
                r = wanted.get(idx)
                if r is None:
                    continue
                dst = os.path.join(dest_dir, r.name + ".pdbqt")
                with open(dst, "w", newline="\n") as out:
                    out.write("\n".join(lines) + "\n")
                items[id(r)] = LigandItem(r.name, dst, src, r.problems)
    return [items[id(r)] for r in records if id(r) in items]


def stage_ligands(inputs, dest_dir, recursive=False, progress=None):
    """Scan *inputs* and stage every ligand found into dest_dir (see stage_records)."""
    return stage_records(scan_ligands(inputs, recursive, progress=progress), dest_dir)


def write_ligand_csv(records, path):
    """Save ligand properties (one row per record) as CSV."""
    import csv
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ligand", "heavy_atoms", "torsions", "dim_1_A", "dim_2_A", "dim_3_A", "longest_A",
                    "charge", "acceptors", "donors_polar_H", "folder", "file", "index_in_file", "smiles",
                    "valid", "problems"])
        for r in records:
            p = r.props
            folder, fname = os.path.split(r.source)
            w.writerow([r.name, p["n_atoms"], "" if p["torsions"] is None else p["torsions"],
                        "%.2f" % p["dims"][0], "%.2f" % p["dims"][1], "%.2f" % p["dims"][2],
                        "%.2f" % p["longest"], "%.3f" % p["charge"], p["n_acceptors"], p["n_polar_h"],
                        folder, fname, r.index + 1, p["smiles"], "no" if r.problems else "yes",
                        "; ".join(r.problems)])
