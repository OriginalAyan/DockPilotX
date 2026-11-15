"""Show docking results in the ChimeraX scene.

* A ligand's whole Vina output (all poses) is opened once as one model (#N with sub-models #N.1, #N.2, ...).
  Selecting a pose shows that sub-model and hides the others. At most *max_open* ligand files stay
  open; the least recently used one is closed.
* H-bonds and clashes are calculated per pose, drawn as pseudobonds in two shared models
  ("Lattice H-bonds", "Lattice clashes") and saved to <workdir>/pose_cache so they come back
  when the pose is shown again, also in a later session.
* Receptor residues revealed for the bonds of a pose are hidden again when that pose is no longer shown.
"""

import os
from collections import OrderedDict

from chimerax.core.commands import run

from .analysis_store import InteractionCache

DEFAULT_COLORS = {"hbond": "#00bfff", "hbond_relaxed": "#2f6fd6", "clash": "#ff6a00"}
SLOP_SENTINEL = "#fe12ab"   # marker colour: bonds drawn in it met only the relaxed criteria
HB_GROUP, CL_GROUP = "Lattice H-bonds", "Lattice clashes"


def hex_to_rgba(text):
    t = (text or "").lstrip("#")
    try:
        return [int(t[0:2], 16), int(t[2:4], 16), int(t[4:6], 16), 255]
    except (ValueError, IndexError):
        return [128, 128, 128, 255]


def _alive(m):
    if m is None:
        return False
    try:
        return not getattr(m, "deleted", False)
    except Exception:
        return False


def _spec(model):
    return "#" + ".".join(str(i) for i in model.id)


def _same_file(a, b):
    try:
        return (os.path.basename(a).lower() == os.path.basename(b).lower()
                and os.path.getsize(a) == os.path.getsize(b))
    except OSError:
        return False


class _Entry:
    def __init__(self, parent, poses):
        self.parent = parent
        self.poses = poses      # list of pose models; poses[k-1] is pose k
        self.bonds = {}         # pose number -> [(pseudobond, kind)], kind in hbond/hbond_relaxed/clash


class PoseViewer:
    def __init__(self, session, host, max_open=5, colors=None):
        """*host* provides get_receptor_models(), set_receptor_models(models, path), close_receptor()."""
        self.session = session
        self.host = host
        self.max_open = max_open
        self.colors = dict(DEFAULT_COLORS)
        self.colors.update(colors or {})
        self.workdir = None
        self.receptor_file = None
        self.cache = None
        self._entries = OrderedDict()   # ligand name -> _Entry, least recently used first
        self._groups = {}               # "hb"/"cl" -> shared pseudobond group
        self._revealed = {}             # receptor residue -> {"owners": set, "hidden": [atoms]}
        self._rec_key = None
        self._rec_atoms = []

    # ---- project ---------------------------------------------------------
    def set_project(self, workdir, receptor_file):
        if workdir != self.workdir:
            self.close_all()
        self.workdir = workdir
        self.receptor_file = receptor_file
        self.cache = InteractionCache(workdir)

    def close_all(self):
        for name in list(self._entries):
            self._evict(name)

    # ---- receptor --------------------------------------------------------
    def _top_models(self):
        return [m for m in self.session.models.list() if len(m.id) == 1]

    def _viewdock_tools(self):
        try:
            return [t for t in self.session.tools.list()
                    if type(t).__name__ == "ViewDockTool" or getattr(t, "tool_name", "") == "ViewDock"]
        except Exception:
            return []

    def _close_new_viewdock(self, before):
        for t in self._viewdock_tools():
            if t not in before:
                try:
                    t.delete()
                except Exception:
                    pass

    def _open(self, path, quiet_viewdock=False):
        """Open a file; with *quiet_viewdock* the ViewDock window that ChimeraX starts for it is closed again."""
        before = set(self._top_models())
        vd_before = list(self._viewdock_tools()) if quiet_viewdock else []
        run(self.session, 'open "%s" format pdbqt' % os.path.abspath(path).replace("\\", "/"))
        if quiet_viewdock:
            self._close_new_viewdock(vd_before)
            try:  # ChimeraX may start ViewDock a moment after the open command returns
                from Qt.QtCore import QTimer
                QTimer.singleShot(300, lambda: self._close_new_viewdock(vd_before))
            except Exception:
                pass
        return [m for m in self._top_models() if m not in before]

    def ensure_receptor(self):
        """Exactly one live copy of the receptor in the scene; never opened twice."""
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

    def _receptor_atoms(self, models):
        key = tuple(id(m) for m in models)
        if key != self._rec_key:
            atoms = []
            for m in models:
                atoms.extend(list(m.atoms))
            self._rec_key, self._rec_atoms = key, atoms
        return self._rec_atoms

    # ---- ligand files ----------------------------------------------------
    def _ligand(self, job):
        """Open the whole output file of *job* once; returns its _Entry."""
        entry = self._entries.get(job.name)
        if entry is not None and _alive(entry.parent):
            self._entries.move_to_end(job.name)
            return entry
        if entry is not None:
            self._evict(job.name)
        new = self._open(job.out_path, quiet_viewdock=True)
        if not new:
            raise RuntimeError("ChimeraX did not open %s" % job.out_path)
        if len(new) == 1:
            parent = new[0]
            kids = sorted(parent.child_models(), key=lambda m: m.id)
            poses = kids if kids else [parent]
        else:
            parent, poses = new[0], sorted(new, key=lambda m: m.id)
        entry = _Entry(parent, poses)
        self._entries[job.name] = entry
        return entry

    def _evict(self, name):
        entry = self._entries.pop(name, None)
        if entry is None:
            return
        for pose in list(entry.bonds):
            self._release((name, pose))
        entry.bonds.clear()
        models = [m for m in {id(m): m for m in ([entry.parent] + entry.poses)}.values() if _alive(m)]
        top = [m for m in models if len(m.id) == 1] or models[:1]
        try:
            self.session.models.close(top)
        except Exception:
            pass

    def pose_model(self, job, pose_no):
        entry = self._ligand(job)
        return entry.poses[min(max(pose_no, 1), len(entry.poses)) - 1]

    def ligand_spec(self, job):
        return _spec(self._ligand(job).parent)

    # ---- showing ---------------------------------------------------------
    def show(self, wanted, with_receptor=True):
        """Show exactly the poses in *wanted* [(job, pose_no)]; hide everything else. Returns shown keys."""
        by_ligand = OrderedDict()
        for job, pose in wanted:
            if job.name not in by_ligand and len(by_ligand) >= max(1, self.max_open):
                continue
            by_ligand.setdefault(job.name, (job, set()))[1].add(pose)
        if with_receptor:
            self.ensure_receptor()
        for job, poses in by_ligand.values():
            self._ligand(job)
        for name, entry in list(self._entries.items()):
            if not _alive(entry.parent):
                self._evict(name)
                continue
            shown = by_ligand[name][1] if name in by_ligand else set()
            entry.parent.display = bool(shown)
            for k, m in enumerate(entry.poses, start=1):
                if _alive(m):
                    m.display = k in shown if m is not entry.parent else bool(shown)
            for k in range(1, len(entry.poses) + 1):
                visible = k in shown
                if visible:
                    self._restore_saved(by_ligand[name][0], k, entry)
                for pb, _kind in entry.bonds.get(k, []):
                    self._set_pb_display(pb, visible)
                if visible:
                    self._reveal((name, k), entry)
                else:
                    self._release((name, k))
        keys = []
        for name, (job, poses) in by_ligand.items():
            self._entries.move_to_end(name)
            keys += [(name, p) for p in sorted(poses)]
        while len(self._entries) > self.max_open:
            victim = next((n for n in self._entries if n not in by_ligand), None)
            if victim is None:
                break
            self._evict(victim)
        return keys

    def hide_all(self):
        self.show([])

    @staticmethod
    def _set_pb_display(pb, visible):
        try:
            pb.display = visible
        except Exception:
            pass

    # ---- revealed receptor residues -----------------------------------------
    def _reveal(self, owner, entry):
        for pb, _kind in entry.bonds.get(owner[1], []):
            try:
                atoms = list(pb.atoms)
            except Exception:
                continue
            for a in atoms:
                try:
                    if a.structure in entry.poses or a.structure is entry.parent:
                        continue
                    r = a.residue
                    st = self._revealed.get(r)
                    if st is None:
                        hidden = [x for x in r.atoms if not x.display]
                        st = self._revealed[r] = {"owners": set(), "hidden": hidden}
                        for x in hidden:
                            x.display = True
                    st["owners"].add(owner)
                except Exception:
                    continue

    def _release(self, owner):
        """Hide again receptor atoms that were only revealed for *owner*'s bonds."""
        for r in list(self._revealed):
            st = self._revealed[r]
            if owner in st["owners"]:
                st["owners"].discard(owner)
                if not st["owners"]:
                    for x in st["hidden"]:
                        try:
                            x.display = False
                        except Exception:
                            pass
                    del self._revealed[r]

    # ---- drawing / colours -------------------------------------------------
    def _group(self, kind):
        g = self._groups.get(kind)
        if not _alive(g):
            name = HB_GROUP if kind == "hb" else CL_GROUP
            g = self.session.pb_manager.get_group(name)
            try:
                if getattr(g, "id", None) is None:
                    self.session.models.add([g])
                g.display = True
                g.dashes = 6
            except Exception:
                pass
            self._groups[kind] = g
        return g

    def _color_of(self, kind):
        return hex_to_rgba(self.colors[kind])

    def set_colors(self, colors):
        self.colors.update(colors)
        for entry in self._entries.values():
            for bonds in entry.bonds.values():
                for pb, kind in bonds:
                    try:
                        pb.color = self._color_of(kind)
                    except Exception:
                        pass

    def _draw(self, entry, pose_no, items, group_kind):
        """Create pseudobonds for *items* [(ligand_atom, receptor_atom, kind)]; replaces that group's bonds for the pose."""
        keep, drop = [], []
        for pb, kind in entry.bonds.get(pose_no, []):
            (drop if (kind == "clash") == (group_kind == "cl") else keep).append((pb, kind))
        for pb, _ in drop:
            try:
                pb.delete()
            except Exception:
                pass
        grp = self._group(group_kind)
        for a, b, kind in items:
            pb = grp.new_pseudobond(a, b)
            try:
                pb.color = self._color_of(kind)
                pb.radius = 0.08
                pb.display = True
            except Exception:
                pass
            keep.append((pb, kind))
        entry.bonds[pose_no] = keep

    def _restore_saved(self, job, pose_no, entry):
        """Redraw saved H-bonds/clashes for a pose that has none drawn yet."""
        if entry.bonds.get(pose_no) or self.cache is None:
            return
        rec = self.cache.load(job.name, pose_no)
        if not rec or not (rec.get("hbonds") or rec.get("clashes")):
            return
        receptors = self.ensure_receptor()
        if not receptors:
            return
        rat = self._receptor_atoms(receptors)
        lat = list(entry.poses[pose_no - 1].atoms)
        if rec.get("n_receptor_atoms") not in (None, len(rat)) or rec.get("n_ligand_atoms") not in (None, len(lat)):
            return  # different structure than when saved
        try:
            hb = [(lat[i], rat[j], "hbond_relaxed" if rel else "hbond") for i, j, rel, _d in rec.get("hbonds", [])]
            cl = [(lat[i], rat[j], "clash") for i, j, _d in rec.get("clashes", [])]
        except IndexError:
            return
        if hb:
            self._draw(entry, pose_no, hb, "hb")
        if cl:
            self._draw(entry, pose_no, cl, "cl")

    # ---- interaction analysis --------------------------------------------
    def _pairs(self, cmd, name, lig_structure, receptors):
        run(self.session, cmd)
        group = self._find_group(name)
        out = []
        if group is not None:
            rec_set = set(receptors)
            try:
                for pb in group.pseudobonds:
                    a, b = pb.atoms
                    try:
                        col = tuple(int(x) for x in pb.color[:3])
                    except Exception:
                        col = None
                    if a.structure is lig_structure and b.structure in rec_set:
                        out.append((a, b, float(pb.length), col))
                    elif b.structure is lig_structure and a.structure in rec_set:
                        out.append((b, a, float(pb.length), col))
            except Exception:
                pass
            try:
                self.session.models.close([group])
            except Exception:
                pass
        return out

    def analyze(self, job, pose_no, kind):
        """Calculate 'hbonds' or 'clashes' between one pose and the receptor, draw and save them.

        Returns (count, relaxed_count, partners): for clashes relaxed_count is 0.
        """
        receptors = self.ensure_receptor()
        if not receptors:
            raise RuntimeError("The receptor could not be opened.")
        entry = self._ligand(job)
        pose = entry.poses[pose_no - 1]
        pose.display = True
        entry.parent.display = True
        lig, rec = _spec(pose), ",".join(_spec(m) for m in receptors)
        tmp = "vpx_tmp_%s" % kind
        rat = self._receptor_atoms(receptors)
        lat = list(pose.atoms)
        li = {a: i for i, a in enumerate(lat)}
        ri = {a: i for i, a in enumerate(rat)}
        if kind == "hbonds":
            # one run: ChimeraX colours bonds that only meet the relaxed criteria with slopColor
            cmd = "hbonds %s,%s restrict %s interModel true intraMol false reveal false select false " \
                  "log false makePseudobonds true twoColors true slopColor %s name %s" % (lig, rec, lig, SLOP_SENTINEL, tmp)
            found = self._pairs(cmd, tmp, pose, receptors)
            slop = tuple(hex_to_rgba(SLOP_SENTINEL)[:3])
            items = [(a, b, "hbond_relaxed" if col == slop else "hbond") for a, b, _, col in found]
            dist = {(a, b): d for a, b, d, _ in found}
            group_kind = "hb"
            saved = {"hbonds": [[li[a], ri[b], 1 if k == "hbond_relaxed" else 0, round(dist.get((a, b), 0.0), 3)]
                                for a, b, k in items]}
        else:
            cmd = "clashes %s restrict %s interModel true intraMol false reveal false select false " \
                  "log false makePseudobonds true name %s" % (lig, rec, tmp)
            found = self._pairs(cmd, tmp, pose, receptors)
            items = [(a, b, "clash") for a, b, _, _c in found]
            group_kind = "cl"
            saved = {"clashes": [[li[a], ri[b], round(d, 3)] for a, b, d, _c in found]}
        self._draw(entry, pose_no, items, group_kind)
        saved.update(n_ligand_atoms=len(lat), n_receptor_atoms=len(rat))
        if self.cache is not None:
            self.cache.save(job.name, pose_no, **saved)
        self._release((job.name, pose_no))
        self._reveal((job.name, pose_no), entry)
        partners = []
        if kind == "hbonds":
            for a, b, k in items:
                partners.append("%s - %s%s" % (self._atom_label(a), self._atom_label(b),
                                               " (relaxed)" if k == "hbond_relaxed" else ""))
        else:
            partners = ["%s - %s" % (self._atom_label(a), self._atom_label(b)) for a, b, _ in items]
        n_relaxed = sum(1 for _, _, k in items if k == "hbond_relaxed")
        return len(items) - n_relaxed, n_relaxed, partners

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
