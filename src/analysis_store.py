"""Persistent per-pose analysis (ratings, H-bond and clash counts) stored in the working folder.

Pure Python. File: <workdir>/analysis.json
"""

import json
import os

from . import layout

FILENAME = "analysis.json"


class AnalysisStore:
    def __init__(self, workdir):
        self.workdir = workdir
        self.path = os.path.join(workdir, FILENAME)
        self.data = {"version": 1, "ligands": {}}
        self.load()

    def load(self):
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                d = json.load(f)
            if isinstance(d, dict) and isinstance(d.get("ligands"), dict):
                self.data = d
        except (OSError, ValueError):
            pass

    def save(self):
        try:
            os.makedirs(self.workdir, exist_ok=True)
            tmp = self.path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(self.data, f, indent=1)
            os.replace(tmp, self.path)
            return True
        except OSError:
            return False

    def get(self, ligand, pose):
        return self.data["ligands"].get(ligand, {}).get(str(pose), {})

    def set(self, ligand, pose, **fields):
        rec = self.data["ligands"].setdefault(ligand, {}).setdefault(str(pose), {})
        for k, v in fields.items():
            if v is None or v == 0 and k == "rating":
                rec.pop(k, None)
            else:
                rec[k] = v
        if not rec:
            self.data["ligands"][ligand].pop(str(pose), None)
        self.save()

    def rated(self, minimum=1):
        """[(ligand, pose_number, record)] with rating >= minimum."""
        out = []
        for lig, poses in self.data["ligands"].items():
            for n, rec in poses.items():
                if int(rec.get("rating", 0) or 0) >= minimum:
                    out.append((lig, int(n), rec))
        return out


class InteractionCache:
    """Saved H-bond / clash atom pairs per pose: <workdir>/pose_cache/<ligand>_pose<N>_interactions.json

    Atoms are stored as indexes into the pose's atom list and the receptor's atom list, so the bonds can
    be redrawn later without recalculating. Layout::

        {"version": 1, "n_ligand_atoms": int, "n_receptor_atoms": int,
         "hbonds":  [[ligand_atom, receptor_atom, relaxed(0/1), distance], ...],
         "clashes": [[ligand_atom, receptor_atom, distance], ...]}
    """

    def __init__(self, workdir):
        self.workdir = workdir
        self.dir = layout.existing(workdir, "cache")

    def path(self, ligand, pose):
        return os.path.join(self.dir, "%s_pose%d_interactions.json" % (ligand, pose))

    def load(self, ligand, pose):
        try:
            with open(self.path(ligand, pose), "r", encoding="utf-8") as f:
                d = json.load(f)
            return d if isinstance(d, dict) else {}
        except (OSError, ValueError):
            return {}

    def save(self, ligand, pose, **fields):
        """Merge *fields* into the stored record and write it atomically."""
        d = self.load(ligand, pose)
        d["version"] = 1
        d.update(fields)
        try:
            if self.dir == layout.path(self.workdir, "cache"):
                layout.make(self.workdir, "cache")
            os.makedirs(self.dir, exist_ok=True)
            p = self.path(ligand, pose)
            with open(p + ".tmp", "w", encoding="utf-8") as f:
                json.dump(d, f)
            os.replace(p + ".tmp", p)
            return True
        except OSError:
            return False
