"""The project file (``project.vinapx``): a small JSON manifest in the working directory.

It records how the receptor was prepared (one *revision* per start from a protein, with the checkpoint
file of every stage), the ligand inputs, the search box and the Vina parameters, so a project can be
reopened later and continued where it stopped. Paths inside the working directory are stored relative to it,
so the folder can be moved. Pure Python.
"""

import datetime
import json
import os

from . import layout

STAGES = ("clean", "hydrogens", "minimized", "pdbqt")
STAGE_FILES = {"clean": "01_clean.pdb", "hydrogens": "02_hydrogens.pdb", "minimized": "03_minimized.pdb"}
FORMAT = "vinapx"


def _now():
    return datetime.datetime.now().replace(microsecond=0).isoformat()


class Project:
    def __init__(self, workdir):
        self.workdir = os.path.abspath(workdir)
        self.path = os.path.join(self.workdir, layout.PROJECT_FILE)
        self.data = self._blank()
        self.loaded = False

    @staticmethod
    def _blank():
        return {"format": FORMAT, "version": 1, "tool_version": "", "created": _now(), "updated": _now(),
                "receptor": {"active": None, "revisions": []},
                "ligands": {"inputs": []}, "box": None, "params": {}, "last_run": None}

    # ---- file ---------------------------------------------------------------
    @classmethod
    def exists_in(cls, workdir):
        return os.path.isfile(os.path.join(workdir, layout.PROJECT_FILE))

    @classmethod
    def locate(cls, path):
        """The project file for a path that may be the file itself or its folder; None if there is none."""
        path = os.path.abspath(os.path.expanduser(str(path).strip().strip('"')))
        if os.path.isdir(path):
            path = os.path.join(path, layout.PROJECT_FILE)
        return path if os.path.isfile(path) else None

    @classmethod
    def open(cls, path):
        p = cls.locate(path)
        if p is None:
            raise FileNotFoundError("No %s found in %s" % (layout.PROJECT_FILE, path))
        proj = cls(os.path.dirname(p))
        proj.load()
        return proj

    def load(self):
        try:
            with open(self.path, encoding="utf-8") as f:
                d = json.load(f)
        except (OSError, ValueError):
            return False
        if not isinstance(d, dict) or d.get("format") != FORMAT:
            return False
        blank = self._blank()
        for k, v in blank.items():
            d.setdefault(k, v)
        self.data, self.loaded = d, True
        return True

    def save(self, tool_version=None):
        if tool_version:
            self.data["tool_version"] = tool_version
        self.data["updated"] = _now()
        try:
            os.makedirs(self.workdir, exist_ok=True)
            tmp = self.path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(self.data, f, indent=2)
            os.replace(tmp, self.path)
            return True
        except OSError:
            return False

    # ---- paths ---------------------------------------------------------------
    def rel(self, p):
        if not p:
            return p
        p = os.path.abspath(p)
        try:
            r = os.path.relpath(p, self.workdir)
        except ValueError:          # another drive on Windows
            return p
        return p if r.startswith("..") else r.replace("\\", "/")

    def abs(self, p):
        if not p:
            return p
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(self.workdir, p))

    # ---- receptor revisions --------------------------------------------------------
    @property
    def revisions(self):
        return self.data["receptor"]["revisions"]

    def revision(self, rev_id):
        return next((r for r in self.revisions if r["id"] == rev_id), None)

    def active(self):
        return self.revision(self.data["receptor"]["active"])

    def new_revision(self, name, source, options=None):
        """Start a new receptor revision; it becomes the active one. Its folder is created."""
        n = 1 + max([int(r["id"][1:]) for r in self.revisions] or [0])
        rev = {"id": "r%03d" % n, "name": name, "source": source, "created": _now(), "options": options or {},
               "stages": {}, "hydrogens": {}, "note": "", "approved": False}
        self.revisions.append(rev)
        self.data["receptor"]["active"] = rev["id"]
        os.makedirs(self.rev_dir(rev), exist_ok=True)
        return rev

    def rev_dir(self, rev):
        base = os.path.join(layout.path(self.workdir, "prepared"), "receptor")
        return os.path.join(base, "%s_%s" % (rev["id"], rev["name"]))

    def stage_file(self, rev, stage):
        return os.path.join(self.rev_dir(rev), STAGE_FILES[stage])

    def set_stage(self, rev, stage, path, **info):
        """Record that *stage* of *rev* finished and its file is *path*. Later stages become invalid for
        stages that depend on this one (minimized needs hydrogens, hydrogens needs clean)."""
        rev["stages"][stage] = dict(info, file=self.rel(path), time=_now())
        rev["approved"] = False         # a changed structure or a new conversion must be approved again
        order = list(STAGES)
        for later in order[order.index(stage) + 1:]:
            rev["stages"].pop(later, None)

    def latest_stage(self, rev):
        """(stage, absolute path) of the last pdb checkpoint that exists, or (None, None)."""
        for stage in ("minimized", "hydrogens", "clean"):
            s = rev["stages"].get(stage)
            if s and os.path.isfile(self.abs(s["file"])):
                return stage, self.abs(s["file"])
        return None, None

    def pdbqt(self, rev=None):
        rev = rev or self.active()
        if not rev:
            return None
        s = rev["stages"].get("pdbqt")
        return self.abs(s["file"]) if s and os.path.isfile(self.abs(s["file"])) else None

    def approve(self, rev, flag=True):
        """Mark the PDBQT of *rev* as approved for docking (only approved ones are offered on the Inputs tab)."""
        rev["approved"] = bool(flag) and self.pdbqt(rev) is not None
        return rev["approved"]

    def approved(self):
        """[(revision, absolute pdbqt path)] of the approved receptors whose file still exists."""
        out = []
        for r in self.revisions:
            p = self.pdbqt(r)
            if r.get("approved") and p:
                out.append((r, p))
        return out

    def summary(self, rev):
        marks = []
        for st in STAGES:
            ok = st in rev["stages"] and os.path.exists(self.abs(rev["stages"][st]["file"]))
            marks.append("%s %s" % (st, "✔" if ok else "–"))
        return "%s  %s  [%s]" % (rev["id"], rev["name"], ", ".join(marks))

    # ---- inputs of the docking run ---------------------------------------------------------
    def set_inputs(self, ligand_paths=None, box=None, params=None, last_run=None):
        if ligand_paths is not None:
            self.data["ligands"]["inputs"] = [self.rel(p) for p in ligand_paths]
        if box is not None:
            c, s = box
            self.data["box"] = {"center": [float(x) for x in c], "size": [float(x) for x in s]}
        if params is not None:
            self.data["params"] = dict(params)
        if last_run is not None:
            self.data["last_run"] = last_run

    def ligand_paths(self):
        return [self.abs(p) for p in self.data["ligands"]["inputs"]]
