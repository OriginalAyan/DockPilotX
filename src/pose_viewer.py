"""Show docking poses in the ChimeraX scene one at a time, with a fixed cap on open models.

A pose is extracted from the Vina output file into its own single-model PDBQT, opened on demand and
shown; poses not currently requested are hidden. When more than *max_open* pose models are loaded the
least recently used one is closed (together with its H-bond/clash pseudobonds), so a large screen never
fills the scene or memory. H-bonds and clashes are calculated only for the pose you ask for.
"""

import os
from collections import OrderedDict

from chimerax.core.commands import run

from .docking import extract_pose


def _spec(model):
    return "#" + ".".join(str(i) for i in model.id)


def _alive(m):
    return m is not None and not getattr(m, "deleted", False)


def _same_file(a, b):
    try:
        return (os.path.basename(a).lower() == os.path.basename(b).lower()
                and os.path.getsize(a) == os.path.getsize(b))
    except OSError:
        return False


class _Entry:
    def __init__(self, model):
        self.model = model
        self.groups = {}  # "hbonds"/"clashes" -> pseudobond group model


class PoseViewer:
    def __init__(self, session, host, max_open=5):
        """*host* provides get_receptor_models(), set_receptor_models(models, path), close_receptor()."""
        self.session = session
        self.host = host
        self.max_open = max_open
        self.workdir = None
        self.receptor_file = None
        self._entries = OrderedDict()   # (ligand, pose) -> _Entry, least recently used first

    # ---- project ---------------------------------------------------------
    def set_project(self, workdir, receptor_file):
        if workdir != self.workdir:
            self.close_all()
        self.workdir = workdir
        self.receptor_file = receptor_file

    def close_all(self):
        for key in list(self._entries):
            self._evict(key)

    # ---- helpers ---------------------------------------------------------
    def _top_models(self):
        return [m for m in self.session.models.list() if len(m.id) == 1]

    def _open(self, path):
        before = set(self._top_models())
        run(self.session, 'open "%s" format pdbqt' % os.path.abspath(path).replace("\\", "/"))
        return [m for m in self._top_models() if m not in before]

    def ensure_receptor(self):
        """Make sure exactly one live copy of the receptor is in the scene; never opens it twice."""
        models = [m for m in self.host.get_receptor_models() if _alive(m)]
        if models:
            current = getattr(self.host, "receptor_path_hint", None)
            if current is None or not self.receptor_file or _same_file(current, self.receptor_file):
                for m in models:
                    m.display = True
                return models
            self.host.close_receptor()
        if not self.receptor_file or not os.path.isfile(self.receptor_file):
            return []
        new = self._open(self.receptor_file)
        self.host.set_receptor_models(new, self.receptor_file)
        return new

    def _load(self, job, pose_no):
        key = (job.name, pose_no)
        entry = self._entries.get(key)
        if entry is not None and _alive(entry.model):
            self._entries.move_to_end(key)
            return entry
        if entry is not None:
            self._evict(key)
        cache = os.path.join(self.workdir, "pose_cache", "%s_pose%d.pdbqt" % (job.name, pose_no))
        if not os.path.isfile(cache):
            extract_pose(job.out_path, pose_no, cache)
        new = self._open(cache)
        if not new:
            raise RuntimeError("ChimeraX did not open %s" % cache)
        model = new[0]
        try:
            model.name = "%s pose %d (%.2f)" % (job.name, pose_no, job.poses[pose_no - 1][0])
        except Exception:
            pass
        entry = _Entry(model)
        self._entries[key] = entry
        return entry

    def _evict(self, key):
        entry = self._entries.pop(key, None)
        if entry is None:
            return
        to_close = [g for g in entry.groups.values() if _alive(g)]
        if _alive(entry.model):
            to_close.append(entry.model)
        if to_close:
            try:
                self.session.models.close(to_close)
            except Exception:
                pass

    # ---- showing ---------------------------------------------------------
    def show(self, wanted, with_receptor=True):
        """Show exactly the poses in *wanted* [(job, pose_no)]; hide the rest. Returns the keys shown."""
        wanted = list(wanted)[: max(1, self.max_open)]
        keys = []
        if with_receptor:
            self.ensure_receptor()
        for job, pose_no in wanted:
            self._load(job, pose_no)
            keys.append((job.name, pose_no))
        for key, entry in self._entries.items():
            visible = key in keys
            if _alive(entry.model):
                entry.model.display = visible
            for g in entry.groups.values():
                if _alive(g):
                    g.display = visible
        # keep the cap: close least-recently-used poses that are not on screen
        while len(self._entries) > self.max_open:
            victim = next((k for k in self._entries if k not in keys), None)
            if victim is None:
                break
            self._evict(victim)
        for k in keys:
            self._entries.move_to_end(k)
        return keys

    def hide_all(self):
        for entry in self._entries.values():
            if _alive(entry.model):
                entry.model.display = False
            for g in entry.groups.values():
                if _alive(g):
                    g.display = False

    # ---- interaction analysis --------------------------------------------
    def analyze(self, job, pose_no, kind):
        """Calculate H-bonds ('hbonds') or clashes ('clashes') between one pose and the receptor.

        Returns (count, partners[list of text]) and leaves the bonds visible in the scene.
        """
        receptors = self.ensure_receptor()
        if not receptors:
            raise RuntimeError("The receptor could not be opened.")
        entry = self._load(job, pose_no)
        entry.model.display = True
        name = "vpx_%s_%s_%d" % ("hb" if kind == "hbonds" else "cl", job.name, pose_no)
        old = entry.groups.pop(kind, None)
        if _alive(old):
            self.session.models.close([old])
        lig, rec = _spec(entry.model), ",".join(_spec(m) for m in receptors)
        if kind == "hbonds":
            cmd = "hbonds %s restrict %s interModel true intraMol false reveal true name %s log false" % (lig, rec, name)
        else:
            cmd = "clashes %s restrict %s interModel true intraMol false reveal true select false name %s log false" % (
                lig, rec, name)
        run(self.session, cmd)
        group = self._find_group(name)
        count, partners = 0, []
        if group is not None:
            entry.groups[kind] = group
            group.display = True
            try:
                for pb in group.pseudobonds:
                    a, b = pb.atoms
                    partners.append("%s - %s" % (self._atom_label(a), self._atom_label(b)))
            except Exception:
                partners = []
            try:
                count = int(group.num_pseudobonds)
            except Exception:
                count = len(partners)
        return count, partners

    def _find_group(self, name):
        try:
            from chimerax.atomic import PseudobondGroup
            groups = self.session.models.list(type=PseudobondGroup)
        except Exception:
            groups = self.session.models.list()
        for g in groups:
            if g.name == name:
                return g
        return None

    @staticmethod
    def _atom_label(atom):
        try:
            r = atom.residue
            return "%s%s%s:%s" % (r.name, r.number, ("." + r.chain_id) if r.chain_id.strip() else "", atom.name)
        except Exception:
            return str(atom)
