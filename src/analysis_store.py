"""Persistent per-pose analysis (ratings, H-bond and clash counts) stored in the working folder.

Pure Python. File: <workdir>/analysis.json
"""

import json
import os

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
