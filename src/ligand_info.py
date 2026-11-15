"""Per-ligand information derived from PDBQT files (for the ligand table).

Pure Python + NumPy; no ChimeraX imports.

Dimensions are computed from atom *centres* (hydrogens kept as written in the
PDBQT), so the true van-der-Waals extent is roughly 3-4 Å larger. They are
orientation independent: the three extents are measured along the ligand's own
principal axes, and the "max span" is the longest atom-to-atom distance.
"""

import os
import re

import numpy as np

from .pdbqt_tools import iter_ligand_blocks, _rec, _atom_type

ELEMENT_OF_TYPE = {
    "C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "S": "S", "SA": "S",
    "H": "H", "HD": "H", "HS": "H", "F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P",
    "MG": "Mg", "MN": "Mn", "FE": "Fe", "ZN": "Zn", "CA": "Ca",
}
ACCEPTOR_TYPES = {"OA", "NA", "SA"}
_SMILES_RE = re.compile(r"^REMARK\s+SMILES\s+(\S+)\s*$")


class LigandInfo:
    """Everything the ligand table shows about one ligand."""
    __slots__ = ("name", "source", "folder", "index", "n_in_file", "n_atoms", "n_heavy", "n_polar_h",
                 "torsions", "hbd", "hba", "charge", "dims", "max_span", "smiles", "problems",
                 "elements", "is_largest")

    def __init__(self):
        self.name = ""
        self.source = ""
        self.folder = ""
        self.index = 1
        self.n_in_file = 1
        self.n_atoms = self.n_heavy = self.n_polar_h = 0
        self.torsions = None
        self.hbd = self.hba = 0
        self.charge = 0.0
        self.dims = (0.0, 0.0, 0.0)
        self.max_span = 0.0
        self.smiles = ""
        self.problems = []
        self.elements = {}
        self.is_largest = False

    @property
    def dims_text(self):
        if self.n_atoms < 2:
            return "-"
        return "%.1f × %.1f × %.1f" % self.dims

    @property
    def formula_text(self):
        order = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
        parts = []
        for el in order + sorted(set(self.elements) - set(order) - {"H"}):
            n = self.elements.get(el)
            if n:
                parts.append(el + (str(n) if n > 1 else ""))
        return "".join(parts)


def _xyz(line):
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except ValueError:
        parts = line.split()
        # fall back: x y z are the three floats before occupancy/bfactor/charge/type
        nums = []
        for p in parts[5:]:
            try:
                nums.append(float(p))
            except ValueError:
                pass
        return nums[0], nums[1], nums[2]


def _charge(line):
    try:
        return float(line[70:76])
    except ValueError:
        return 0.0


def extents_and_span(coords):
    """Return ((L, W, H) along principal axes sorted descending, longest atom-atom distance)."""
    xyz = np.asarray(coords, dtype=float)
    n = len(xyz)
    if n < 2:
        return (0.0, 0.0, 0.0), 0.0
    centred = xyz - xyz.mean(axis=0)
    if n >= 3:
        _, _, vt = np.linalg.svd(centred, full_matrices=False)
        proj = centred @ vt.T
        ext = proj.max(axis=0) - proj.min(axis=0)
        if len(ext) < 3:
            ext = np.concatenate([ext, np.zeros(3 - len(ext))])
    else:
        ext = np.array([np.linalg.norm(xyz[1] - xyz[0]), 0.0, 0.0])
    ext = tuple(float(v) for v in sorted(ext, reverse=True))
    if n <= 1500:
        d = xyz[:, None, :] - xyz[None, :, :]
        span = float(np.sqrt((d * d).sum(-1)).max())
    else:
        span = ext[0]
    return ext, span


def analyze_lines(lines, name="", source=""):
    info = LigandInfo()
    info.name = name
    info.source = source
    info.folder = os.path.dirname(source) if source else ""
    coords = []
    for line in lines:
        rec = _rec(line)
        if rec in ("ATOM", "HETATM"):
            t = _atom_type(line).upper()
            el = ELEMENT_OF_TYPE.get(t, t.title() if t else "?")
            info.n_atoms += 1
            info.elements[el] = info.elements.get(el, 0) + 1
            if el == "H":
                info.n_polar_h += 1
                if t == "HD":
                    info.hbd += 1
            else:
                info.n_heavy += 1
            if t in ACCEPTOR_TYPES:
                info.hba += 1
            info.charge += _charge(line)
            try:
                coords.append(_xyz(line))
            except (IndexError, ValueError):
                pass
        elif rec == "TORSDOF":
            parts = line.split()
            if len(parts) > 1 and parts[1].isdigit():
                info.torsions = int(parts[1])
        elif rec == "REMARK" and not info.smiles:
            m = _SMILES_RE.match(line.strip())
            if m:
                info.smiles = m.group(1)
    if coords:
        info.dims, info.max_span = extents_and_span(coords)
    return info


def scan_ligands(inputs, recursive=False, progress=None, cancelled=None):
    """Yield a LigandInfo for every ligand in files/folders/libraries (streaming).

    progress(file_index, n_files) is called as files are passed.
    """
    last = -1
    for blk in iter_ligand_blocks(inputs, recursive):
        if cancelled is not None and cancelled():
            return
        info = analyze_lines(blk.lines, blk.name, blk.source)
        info.index = blk.index
        info.n_in_file = blk.n_in_file
        info.problems = list(blk.problems)
        if progress and blk.file_index != last:
            last = blk.file_index
            progress(blk.file_index + 1, blk.n_files)
        yield info


def mark_largest(infos):
    """Flag the ligand with the greatest max span; returns it (or None)."""
    best = None
    for i in infos:
        i.is_largest = False
        if not i.problems and (best is None or i.max_span > best.max_span):
            best = i
    if best is not None and best.max_span > 0:
        best.is_largest = True
        return best
    return None


def suggested_box_edge(max_span, padding=4.0):
    """Cubic box edge able to hold the ligand rotated any way, plus padding each side."""
    return float(max_span) + 2.0 * float(padding)
