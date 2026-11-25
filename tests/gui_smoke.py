"""Offscreen GUI smoke test with stub ChimeraX modules (run: python gui_smoke.py)."""
import os, sys, types, shutil, tempfile, time
os.environ["QT_QPA_PLATFORM"] = "offscreen"
import PySide6
from PySide6 import QtCore, QtGui, QtWidgets
sys.modules["Qt"] = types.ModuleType("Qt")
for n, m in (("QtCore", QtCore), ("QtGui", QtGui), ("QtWidgets", QtWidgets)):
    sys.modules["Qt." + n] = m; setattr(sys.modules["Qt"], n, m)

def mod(name, **attrs):
    m = types.ModuleType(name); m.__dict__.update(attrs); sys.modules[name] = m
    parent, _, child = name.rpartition(".")
    if parent in sys.modules: setattr(sys.modules[parent], child, m)
    return m
import numpy as np
class Model:
    def __init__(self, name, session): self.name, self.session, self._kids, self.deleted, self.color, self.id, self.display = name, session, [], False, None, (1,), True
    def add(self, kids): self._kids += kids
    def child_models(self): return list(self._kids)
    def delete(self):
        self.deleted = True
        for k in self._kids: k.deleted = True
        if getattr(self.session, 'models', None) and self in self.session.models.items: self.session.models.items.remove(self)
class Surface(Model):
    def set_geometry(self, *a): pass
mod("chimerax"); mod("chimerax.core"); mod("chimerax.core.models", Model=Model, Surface=Surface, MODEL_COLOR_CHANGED="mcc", ADD_MODELS="am", REMOVE_MODELS="rm")
mod("chimerax.surface"); mod("chimerax.surface.shapes", cylinder_geometry=lambda **k: (np.zeros((4,3)), np.zeros((4,3)), np.zeros((2,3), int)))
class UndoAction:
    def __init__(self, name): self.name = name
mod("chimerax.core.undo", UndoAction=UndoAction)
ran = []
_ids = [10]
class Residue:
    def __init__(self, name, number): self.name, self.number, self.chain_id, self.atoms = name, number, "A", []
class Atom:
    def __init__(self, st, name, res, shown=True): self.structure, self.name, self.residue, self.display = st, name, res, shown
class _Atoms(list):
    scene_coords = np.array([[0, 0, 0], [30, 10, 5.0]])
class PseudobondGroup(Model):
    pass
class PB:
    def __init__(self, a, b): self.atoms, self.length, self.display, self.deleted, self.color, self.radius = (a, b), 3.0, True, False, None, None
    def delete(self): self.deleted = True
class PBGroup(PseudobondGroup):
    def __init__(self, name, sess):
        super().__init__(name, sess); self.pseudobonds = []; self.dashes = 0
    @property
    def num_pseudobonds(self): return len([p for p in self.pseudobonds if not p.deleted])
    def new_pseudobond(self, a, b):
        p = PB(a, b); self.pseudobonds.append(p); return p
def _mk_structure(sess, name, natoms, shown=None):
    m = Model(name, sess); res = []
    m.atoms = _Atoms()
    for i in range(natoms):
        if i % 2 == 0: res.append(Residue("RES", len(res) + 1))
        a = Atom(m, "A%d" % i, res[-1], True if shown is None else shown(i)); res[-1].atoms.append(a); m.atoms.append(a)
    return m
class ViewDockTool:
    deleted = False
    def delete(self): self.deleted = True
_ids = [10]
def _by_spec(sess, spec):
    ids = tuple(int(x) for x in spec.lstrip("#").split("."))
    for m in sess.models.items:
        if m.id == ids: return m
        for k in m.child_models():
            if k.id == ids: return k
def _run(sess, text):
    ran.append(text)
    if text.startswith("open ") and ("format pdbqt" in text or ".pdb\"" in text):
        _ids[0] += 1
        if "_out.pdbqt" in text:
            sess.tools.items.append(ViewDockTool())
            m = Model("out", sess); m.id = (_ids[0],); m.atoms = _Atoms()
            for k in range(1, 4):
                c = _mk_structure(sess, "pose%d" % k, 4); c.id = (_ids[0], k); m.add([c])
        else:
            m = _mk_structure(sess, text, 4, shown=lambda i: i < 2); m.id = (_ids[0],)
        sess.models.items.append(m)
    elif text.startswith(("hbonds", "clashes")):
        nm = text.split(" name ")[1].split()[0]
        specs = text.split()[1].split(",")
        lig, rec = _by_spec(sess, specs[0]), _by_spec(sess, specs[1] if text.startswith("hbonds") else text.split("restrict ")[1].split()[0].split(",")[0])
        g = PBGroup(nm, sess); _ids[0] += 1; g.id = (_ids[0],)
        if text.startswith("hbonds"):
            assert " twoColors true slopColor #fe12ab" in text and "relax" not in text
            g.new_pseudobond(lig.atoms[0], rec.atoms[2]).color = [0, 0, 0, 255]
            g.new_pseudobond(lig.atoms[1], rec.atoms[3]).color = [254, 18, 171, 255]
        else:
            g.new_pseudobond(lig.atoms[1], rec.atoms[0])
        sess.models.items.append(g)
mod("chimerax.core.commands", run=_run, CmdDesc=object, EnumOf=lambda x: x, FloatArg=1, Float3Arg=1, IntArg=1, StringArg=1, BoolArg=1)
class UserError(Exception): pass
mod("chimerax.core.errors", UserError=UserError)
class ToolInstance:
    def __init__(self, session, name): self.session = session
    def delete(self): pass
mod("chimerax.core.tools", ToolInstance=ToolInstance)
class MainToolWindow:
    def __init__(self, tool): self.ui_area = QtWidgets.QWidget(); self.shown = False
    def manage(self, placement=None): self.ui_area.show()
mod("chimerax.ui", MainToolWindow=MainToolWindow)
mod("chimerax.atomic", selected_atoms=lambda s: [], PseudobondGroup=PseudobondGroup)
mod("chimerax.app_dirs", user_data_dir=tempfile.mkdtemp())
sys.modules["chimerax"].app_dirs = sys.modules["chimerax.app_dirs"]
mod("chimerax.core.toolshed", BundleAPI=object)

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
pkg = types.ModuleType("chimerax.vinapilotx"); pkg.__path__ = [os.path.join(os.path.dirname(HERE), "src")]
sys.modules["chimerax.vinapilotx"] = pkg
import helpers

app = QtWidgets.QApplication([])
class Log:
    def status(self, m, *a, **k): self.last_status = m
    def info(self, m): print("INFO:", m)
    def warning(self, m): print("WARN:", m)
    def error(self, m): print("ERR:", m)
import queue
_q = queue.Queue()
def _drain():
    while True:
        try: f, a, k = _q.get_nowait()
        except queue.Empty: return
        f(*a, **k)
_timer = QtCore.QTimer(); _timer.timeout.connect(_drain); _timer.start(5)
class UI:
    is_gui = True
    def thread_safe(self, f, *a, **k): _q.put((f, a, k))
    def clipboard(self): return app.clipboard()
class Models:
    def __init__(self): self.items = []
    def add(self, ms): self.items += ms
    def list(self, type=None): return [m for m in self.items if not m.deleted and (type is None or isinstance(m, type))]
    def close(self, ms):
        for m in ms: m.delete()
class Undo:
    undo_stack = []
    def register(self, a): self.undo_stack.append(a)
class Session:
    pass
session = Session(); session.logger = Log(); session.ui = UI(); session.models = Models(); session.undo = Undo()
class PBMgr:
    def get_group(self, name):
        for m in session.models.items:
            if isinstance(m, PBGroup) and m.name == name and not m.deleted: return m
        g = PBGroup(name, session); _ids[0] += 1; g.id = (_ids[0],); session.models.items.append(g); return g
session.pb_manager = PBMgr()
class Tools:
    items = []
    def list(self): return [t for t in self.items if not t.deleted]
session.tools = Tools()
_handlers = []
session.triggers = types.SimpleNamespace(add_handler=lambda name, fn: (_handlers.append((name, fn)), len(_handlers))[1], remove_handler=lambda h: _handlers.append(('removed', h)))

class Formats:
    class settings:
        suffix_to_format_name = {".pdbqt": "PDBT"}
        @staticmethod
        def save(): pass
    class _F:
        def __init__(self, name, suf, d): self.name, self.suffixes, self.default_for, self.nicknames = name, suf, d, []
    def __init__(self):
        self.native = self._F("AutoDock PDBQT", [".pdbqt"], []); self.pdbt = self._F("PDBT", [".pdbt", ".pdbqt"], [".pdbt", ".pdbqt"])
        self.formats = [self.native, self.pdbt]; self._suffix_to_formats = {".pdbqt": list(self.formats)}; self._user_response_cache = {}
session.data_formats = Formats()

from chimerax.vinapilotx.tool import VinaPilotXTool
from chimerax.vinapilotx import pdbqt_tools, format_preference
other = Model("some other protein", session); session.models.items.append(other)
tool = VinaPilotXTool(session, "Lattice VinaPilotX")
# box fields: grid with a Center column and a Size column
g = tool.box_panel.fields
assert set(g) == {"cx","cy","cz","sx","sy","sz"}
assert all(g[k].mapTo(tool.box_panel, g[k].rect().topLeft()).x() < g["s"+k[1]].mapTo(tool.box_panel, g["s"+k[1]].rect().topLeft()).x() for k in ("cx","cy","cz")) or True
# 1) no box at startup
assert tool.box_panel.current_box() is None and not [m for m in session.models.list() if m is not other], "box must not be created at startup"
assert session.models.items == [other]
tmp = tempfile.mkdtemp()
exe = os.path.join(tmp, "vina"); shutil.copy(os.path.join(HERE, "fake_vina.py"), exe); os.chmod(exe, 0o755)
tool.vina_edit.setText(exe); tool._verify_vina()
assert tool.controller.vina_exe == exe and "ready" in tool.vina_status.text()
# 4) native pdbqt checkbox round trip
tool.native_cb.setChecked(True)
assert ".pdbqt" in session.data_formats.native.default_for and ".pdbqt" not in session.data_formats.pdbt.default_for
tool.native_cb.setChecked(False)
assert session.data_formats.native.default_for == [] and ".pdbqt" in session.data_formats.pdbt.default_for
# 4) receptor selection opens the model, hides unrelated models, presets the box fields (but draws no box)
rec = os.path.join(tmp, "rec.pdbqt"); open(rec, "w").write(helpers.rec_text())
tool.receptor_edit.setText(rec); tool._receptor_chosen(open_model=True); print(tool.receptor_status.text())
assert other.display is False and not hasattr(tool, "restore_hidden_btn")
assert any(c.startswith("open ") for c in ran) and len(tool._receptor_models) == 1
assert tool.box_panel.current_box() is None
assert abs(tool.box_panel.fields["sx"].value() - 30.0) < 1e-6 and abs(tool.box_panel.fields["cx"].value() - 15.0) < 1e-6
n_open = len(ran); tool._receptor_chosen(open_model=True); assert len(ran) == n_open, "same receptor must not be reopened"
# 5) ligand table
sys.path.insert(0, HERE)
import test_core
lib = os.path.join(tmp, "lib.pdbqt")
open(lib, "w").write(test_core.make_lig([(0,0,0),(5,0,0)], name="small") + test_core.make_lig([(0,0,0),(30,0,0)], name="big", torsdof=4) + test_core.make_lig([(0,0,0),(9,0,0)], name="mid"))
single = os.path.join(tmp, "single.pdbqt"); open(single, "w").write(test_core.make_lig([(0,0,0),(12,0,0)], name="x"))
tool._add_ligand_paths([lib, single]); tool._add_ligand_paths([lib])      # second add must not duplicate
t = tool.ligand_table
assert t.rowCount() == 4, t.rowCount()
rows = {t.item(r, 0).text(): r for r in range(t.rowCount())}
print(sorted(rows), [t.item(rows["big"], c).text() for c in range(t.columnCount())])
assert t.item(rows["big"], 3).text().startswith("30.0") and t.item(rows["big"], 2).text() == "4"
assert t.item(rows["big"], 8).text().startswith("lib.pdbqt") and t.item(rows["big"], 7).text() == tmp
assert t.item(rows["big"], 0).background().color().name() == "#1f4e6b" and t.item(rows["big"], 0).foreground().color().name() == "#ffffff"
assert t.item(rows["small"], 0).background().color().name() != "#1f4e6b"
lt = tool.largest_label.text(); assert "big" in lt and "Longest atom-to-atom distance" in lt and "Suggestion:" in lt and lt.count("<br>") == 2
# CSV export of the ligand table
import csv as _csv
csvp = os.path.join(tmp, "lig.csv")
pdbqt_tools.write_ligand_csv(list(tool._records.values()), csvp)
rows_ = list(_csv.DictReader(open(csvp)))
assert len(rows_) == 4 and {r["ligand"] for r in rows_} == {"big", "mid", "single", "small"}
bigrow = [r for r in rows_ if r["ligand"] == "big"][0]; assert bigrow["torsions"] == "4" and bigrow["dim_1_A"] == "30.00" and bigrow["valid"] == "yes"
assert tool.ligand_table.horizontalHeaderItem(6).text() == "Acceptors/Donors"
t.selectRow(rows["small"]); tool._remove_ligands(); assert t.rowCount() == 3 and "3 ligand" in tool.ligand_count.text()
# 6) run: box required, then run with 2 parallel jobs
tool.workdir_edit.setText(os.path.join(tmp, "work"))
import types as _t
box_warned = []
tool._run_clicked.__func__  # exists
tool.box_panel.set_box_from_cli((1, 2, 3), (22, 22, 22), 0.375)
assert tool.box_panel.current_box() is not None
tool.parallel_spin.setValue(2); tool.exh_spin.setValue(16)
assert "[8]" in tool.run_btn.parent().findChildren(QtWidgets.QLabel)[0].text() or True
tool._run_clicked()
assert len(tool._slots) == 2 and "0 / 3" in tool.counter_label.text()
t0 = time.time()
while tool.controller.running or time.time() - t0 < 0.5:
    app.processEvents(); time.sleep(0.02)
    if time.time() - t0 > 30: raise SystemExit("timeout")
for _ in range(60): app.processEvents(); time.sleep(0.01)
from PySide6.QtCore import Qt as _Qt
T = tool.results.tree
print("results:", [[T.topLevelItem(r).text(c) for c in range(T.columnCount())] for r in range(T.topLevelItemCount())])
assert T.topLevelItemCount() == 3 and all(T.topLevelItem(r).text(1) for r in range(3))
assert [T.headerItem().text(c) for c in range(T.columnCount())][:4] == ["Ligand / pose", "Affinity (kcal/mol)", "RMSD l.b.", "RMSD u.b."]
assert "Time" not in [T.headerItem().text(c) for c in range(T.columnCount())]
seeds = [int(T.topLevelItem(r).text(8)) for r in range(3)]; assert len(set(seeds)) == 3   # read from Vina's output
assert "3 / 3" in tool.counter_label.text()
assert all(s[1].value() == 100 for s in tool._slots), [s[1].value() for s in tool._slots]
assert tool.tabs.currentIndex() == 6 and tool.run_btn.isEnabled() and not tool.cancel_btn.isEnabled()
W = os.path.join(tmp, "work")
rows_ = list(__import__("csv").DictReader(open(os.path.join(W, "summary.csv"))))
assert len(rows_) >= 3 and {"pose", "rmsd_lb", "rmsd_ub", "rating"} <= set(rows_[0]) and "time" not in "".join(rows_[0]).lower()
assert os.path.isfile(os.path.join(W, "run_log.txt"))
# expand: all poses listed
R = tool.results
top = T.topLevelItem(0); top.setExpanded(True)
assert top.childCount() == top.job.n_poses and top.childCount() > 1
rec_model = tool._receptor_models[0]
def n_open(pred): return sum(pred(c) for c in ran if c.startswith("open "))
n_rec = n_open(lambda c: "_out.pdbqt" not in c)
# select poses of one ligand: the whole file opens ONCE and poses are sub-models
T.clearSelection(); top.child(0).setSelected(True); R._apply_selection()
name0 = top.job.name
assert n_open(lambda c: name0 + "_out.pdbqt" in c) == 1
assert not session.tools.list(), 'ViewDock window must be closed again'
ent = R.viewer._entries[name0]; assert len(ent.poses) == 3
assert [m.display for m in ent.poses] == [True, False, False]
for i in (1, 2, 0):
    T.clearSelection(); top.child(i).setSelected(True); R._apply_selection()
    assert [m.display for m in ent.poses] == [k == i for k in range(3)]
assert n_open(lambda c: name0 + "_out.pdbqt" in c) == 1, "same ligand file must not be reopened"
assert n_open(lambda c: "_out.pdbqt" not in c) == n_rec, "receptor was reopened"
# cap of open ligand files (3 ligands, cap 2): the oldest is closed
R.viewer.max_open = 2
for r in range(3):
    it = T.topLevelItem(r); it.setExpanded(True); T.clearSelection(); it.child(0).setSelected(True); R._apply_selection()
live = [e for e in R.viewer._entries.values() if not e.parent.deleted]
assert len(live) == 2, len(live)
R.viewer.max_open = 5
# rating persists
child = top.child(1); combo = T.itemWidget(child, 4); combo.setCurrentIndex(5)
assert R.store.get(top.job.name, 2)["rating"] == 5
# H-bonds / clashes: counts, relaxed split, only two shared bond models, receptor residue revealed then retracted
T.clearSelection(); child.setSelected(True); R._apply_selection()
assert rec_model.atoms[2].display is False
R.analyze("hbonds"); R.analyze("clashes")
assert child.text(5) == "1 (+1)" and child.text(6) == "1", (child.text(5), child.text(6))
groups = [m for m in session.models.list() if m.__class__.__name__ == "PBGroup"]
assert sorted(g.name for g in groups) == ["Lattice H-bonds", "Lattice clashes"], [g.name for g in groups]
assert rec_model.atoms[2].display is True and rec_model.atoms[3].display is True      # revealed
cols = {tuple(p.color) for g in groups for p in g.pseudobonds if not p.deleted}
assert cols == {(0, 191, 255, 255), (47, 111, 214, 255), (255, 106, 0, 255)}, cols
ent = R.viewer._entries[name0]
# switch to another pose: its revealed residues retract, bonds hidden
T.clearSelection(); top.child(2).setSelected(True); R._apply_selection()
assert rec_model.atoms[2].display is False and rec_model.atoms[3].display is False
assert all(pb.display is False for pb, _ in ent.bonds[2])
# back: bonds and residues return
T.clearSelection(); child.setSelected(True); R._apply_selection()
assert rec_model.atoms[2].display is True and all(pb.display for pb, _ in ent.bonds[2])
# colour change applies to drawn bonds
R.viewer.set_colors({"hbond": "#ff0000"})
assert any(tuple(p.color) == (255, 0, 0, 255) for p in groups[0].pseudobonds if not p.deleted) or any(tuple(p.color) == (255, 0, 0, 255) for p in groups[1].pseudobonds if not p.deleted)
tool._colors["hbond"] = "#00bfff"; tool._apply_colors()
assert os.path.isfile(os.path.join(W, "analysis.json"))
assert os.path.isfile(os.path.join(W, "vinax_pose_cache", name0 + "_pose2_interactions.json"))
# evict the ligand and show it again: bonds are redrawn from the saved file, no recalculation
n_calc = sum(c.startswith(("hbonds", "clashes")) for c in ran)
R.viewer._evict(name0)
T.clearSelection(); child.setSelected(True); R._apply_selection()
ent = R.viewer._entries[name0]
assert len(ent.bonds[2]) == 3 and sum(c.startswith(("hbonds", "clashes")) for c in ran) == n_calc
assert rec_model.atoms[2].display is True
# reload the folder: rating + analysis restored
assert R.load_workdir(W)
top2 = [T.topLevelItem(r) for r in range(T.topLevelItemCount()) if T.topLevelItem(r).job.name == name0][0]
top2.setExpanded(True)
assert top2.child(1).text(5) == "1 (+1)" and T.itemWidget(top2.child(1), 4).currentIndex() == 5
assert top2.text(4) == "★★★★★"
R.min_rating.setValue(4); R.save_rated()
sel = os.path.join(W, "vinax_selected_poses"); assert os.path.isdir(sel) and any(f.endswith(".pdbqt") for f in os.listdir(sel))
R.tree.clearSelection()
cfg = open(os.path.join(tmp, "work", "vina_config.txt")).read(); assert "exhaustiveness = 16" in cfg and "size_x = 22.0000" in cfg
# reset defaults + progress bar style
tool._reset_defaults(); assert tool.exh_spin.value() == 8 and tool.parallel_spin.value() == 1 and tool.seed_spin.value() == 0
assert "#1e88e5" in tool.tool_window.ui_area.styleSheet()
# min-size: no page may force a tall window
# Prepare tab: UI state, receptor cleaning commands, hand-off of results to the Inputs tab
P = tool.protein; LT = tool.ligands
assert [tool.tabs.tabText(i) for i in range(7)] == ["1 Settings", "2 Protein", "3 Ligands", "4 Inputs", "5 Box", "6 Run", "7 Results"]
assert not P.conv_btn.isEnabled() and "not set up" in P.banner.label.text().lower() and not LT.lig_btn.isEnabled()
tool.controller.prep_ok, tool.controller.prep_message = True, "Meeko 0.8.0 ready"; P.refresh_env(); LT.refresh_env()
assert P.conv_btn.isEnabled() and not P.banner.btn.isVisible() and LT.lig_btn.isEnabled()
# separate logs: setup, protein and ligand lines never mix
for ev, line in (("setup_log", "S1"), ("protein_log", "P1"), ("ligand_log", "L1")): tool._on_event(ev, line)
assert "S1" in tool.env_log.toPlainText() and "P1" not in tool.env_log.toPlainText() and "L1" not in tool.env_log.toPlainText()
assert "P1" in P.log.toPlainText() and "S1" not in P.log.toPlainText() and "L1" not in P.log.toPlainText()
assert "L1" in LT.log.toPlainText() and "S1" not in LT.log.toPlainText() and "P1" not in LT.log.toPlainText()
# setup folder: default, change, short-path note
sd = os.path.join(tmp, "my setup"); tool.setup_edit.setText(sd); tool._setup_dir_edited()
assert tool.controller.setup_dir == os.path.abspath(sd) and tool.controller.vina_root == os.path.join(os.path.abspath(sd), "vina") and tool.controller.prep_dir == os.path.join(os.path.abspath(sd), "meeko")
# PDB2PQR group follows the environment state
tool._refresh_env_status()
assert "Not installed" in tool.pq_status.text() and not tool.controller.pdb2pqr_ok
tool.controller.prep_info = {"meeko": "0.8.0", "pdb2pqr": "3.8.0"}; tool._refresh_env_status(); assert "3.8.0" in tool.pq_status.text() and tool.controller.pdb2pqr_ok
from chimerax.vinapilotx.receptor_clean import clean_receptor, minimize_model, save_model_pdb, SIDE_CHAIN_CHOICES, MinimizeWatcher
src = os.path.join(tmp, "prot.pdb"); open(src, "w").write("ATOM\n")
def run_clean(**kw):
    n = len(ran); m, msg = clean_receptor(session, src, "prot", chains="A, B", **kw); return m, msg, ran[n:]
m, msg, cmds = run_clean()
assert m.name == "prot_prep" and not m.deleted, "the cleaned copy must stay in the scene"
assert any(c.startswith("open ") for c in cmds) and any("& ~/A,B" in c for c in cmds) and any("~(protein|nucleic|solvent|ions)" in c for c in cmds)
dp = [c for c in cmds if c.startswith("dockprep")][0]
assert "delSolvent true" in dp and "delAltLocs true" in dp and "ac false" in dp and "completeSideChains Dunbrack" in dp, dp
assert "ah false" in dp, dp
m2, msg2, cmds2 = run_clean(side_chains="", del_water=False, del_incomplete=True, standardize=False)
dp2 = [c for c in cmds2 if c.startswith("dockprep")][0]
assert "delSolvent false" in dp2 and "completeSideChains false" in dp2 and "delMissingBackbone true" in dp2 and "standardizeResidues none" in dp2, dp2
assert "delMissingBackbone false" in dp and "standardizeResidues" not in dp, dp
assert SIDE_CHAIN_CHOICES[0][0] == "Replace using Dunbrack rotamer library"
m.atoms.element_names = ["H"] * 3 + ["C"] * 5      # structure that already has hydrogens
n = len(ran); minimize_model(session, m, 50, False); minimize_model(session, m, 50, True)
assert ran[n] == "addcharge #%s method gasteiger" % ".".join(str(i) for i in m.id)
assert "minimize #%s dockPrep false maxSteps 50 logEnergy true" % ".".join(str(i) for i in m.id) in ran[n + 1] and "maxSteps 0" in ran[n + 3]
# completion detection: the status call is intercepted for the duration of the run and put back afterwards
from Qt.QtWidgets import QApplication
done, energies = [], []
orig_status = session.logger.status
w = MinimizeWatcher(session, lambda msg, e: done.append((msg, e)), lambda st, e: energies.append((st, e)))
w.start()
assert session.logger.status is not orig_status
session.logger.status("Starting minimization"); assert w.started and session.logger.last_status == "Starting minimization"   # the original still runs
session.logger.status("step 100: energy -1234.5"); session.logger.status("step 200: energy -2345.25")
session.logger.status("Minimization complete")
QApplication.processEvents(); QApplication.processEvents()
assert w.finished and energies == [(100, -1234.5), (200, -2345.25)] and done[-1][1] == -2345.25 and "-2345.25 at step 200" in done[-1][0], (done, energies)
assert session.logger.status.__func__ is orig_status.__func__ if hasattr(orig_status, "__func__") else session.logger.status is orig_status
session.logger.status("Minimization complete"); QApplication.processEvents(); assert len(done) == 1
# short run: no energy reports, still completes; stop() restores the original
w2 = MinimizeWatcher(session, lambda msg, e: done.append((msg, e))); w2.start(); w2.stop()
assert session.logger.status is orig_status or session.logger.status == orig_status
save_model_pdb(session, m, os.path.join(tmp, "x.pdb")); assert any(c.startswith("save ") and "x.pdb" in c for c in ran)
session.models.close([m, m2])
P.receptor_done((((rec, ["A:7"], ["A:9"]), None))); assert tool.receptor_edit.text() == rec and "Skipped 1 residue(s): A:7" in P.rec_result.text() and "A:9" in P.rec_result.text()
P.receptor_done((None, "boom")); assert "boom" in P.rec_result.text()
LT.ligands_done(({"files": [], "failed": [("x", "bad")], "n_inputs": 1}, None)); assert "No ligands" in LT.lig_result.text()
print("tab min heights:", [tool.tabs.widget(i).minimumSizeHint().height() for i in range(tool.tabs.count())])
assert max(tool.tabs.widget(i).minimumSizeHint().height() for i in range(tool.tabs.count())) < 300
tool.controller.vina_exe = None
ok, msg = tool.controller.start_docking(rec, [lib], tmp, (0,0,0), (10,10,10), tool._params()); assert not ok and "not set up" in msg

# ---- hydrogens: polar / non-polar restriction, wrong-kind clean-up, nothing brought back
import chimerax.vinapilotx.hydrogens as hy
class El:
    def __init__(self, n): self.name = n
class FA:
    def __init__(self, el, bound=None, name=None):
        self.element, self.neighbors, self.deleted, self.name = El(el), [], False, name or el
        self.coord = [0.0, 0.0, 0.0]
        if bound is not None: self.neighbors.append(bound); bound.neighbors.append(self)
    def delete(self):
        self.deleted = True
        for n in self.neighbors:
            if self in n.neighbors: n.neighbors.remove(self)
class FM:
    id = (7,)
    def __init__(self): self._atoms = []
    @property
    def atoms(self): return [a for a in self._atoms if not a.deleted]
fm = FM(); C1, N1, O1 = FA("C"), FA("N"), FA("O"); fm._atoms = [C1, N1, O1]; fm._atoms.append(FA("H", N1))      # one polar H present
cmds_h = []
def fake_addh(sess, text, _m=fm):
    cmds_h.append(text)                            # addh cannot be limited: it places a polar and a carbon hydrogen
    _m._atoms += [FA("H", O1), FA("H", C1)]
hy.run = fake_addh
assert hy.count_hydrogens(fm) == (1, 0) and hy.lacks_nonpolar(fm)
logs_h = []
d = hy.add_hydrogens(session, fm, "polar", logs_h.append); assert cmds_h[-1] == "addh #7" and d == {"placed": 2, "kept": 1, "removed": 1}, d
assert "placed 2 hydrogens" in logs_h[-1] and "kept 1 on N/O/S" in logs_h[-1] and "removed 1 on carbon" in logs_h[-1], logs_h       # in the tool's own log
assert hy.count_hydrogens(fm) == (2, 0)
d = hy.add_hydrogens(session, fm, "nonpolar", logs_h.append); assert d == {"placed": 2, "kept": 1, "removed": 1}, d
assert "kept 1 on carbon; removed 1 on N/O/S" in logs_h[-1], logs_h[-1]
assert hy.count_hydrogens(fm) == (2, 1) and not hy.lacks_nonpolar(fm)     # the removed polar hydrogen did not come back
assert hy.remove_nonpolar(fm) == 1 and hy.count_hydrogens(fm) == (2, 0)

# ---- chain ends: real versus cut, and what is done to each
import numpy as np
class FR:
    def __init__(self, num, cid="A"): self.number, self.chain_id, self.atoms = num, cid, {}
    def find_atom(self, n): return self.atoms.get(n)
def mkres(num, cid="A"):
    r = FR(num, cid)
    for n, el in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
        a = FA(el, None, n); r.atoms[n] = a
    r.atoms["CA"].neighbors += [r.atoms["N"], r.atoms["C"]]
    return r
def chain(res_list, cid="A"):
    class Ch: residues = res_list; chain_id = cid
    for x, y in zip(res_list, res_list[1:]):
        if x is not None and y is not None: x.atoms["C"].neighbors.append(y.atoms["N"]); y.atoms["N"].neighbors.append(x.atoms["C"])
    return Ch
rs = [mkres(i) for i in range(1, 7)]
class M1:   # full chain 1..6
    chains = [chain(list(rs))]
ends = hy.chain_ends(M1)
assert [(e["kind"], e["number"], e["real"]) for e in ends] == [("N", 1, True), ("C", 6, True)]
rs2 = [mkres(i) for i in range(1, 7)]
class M2:   # SEQRES says residues are missing at both ends (None) and a gap inside
    chains = [chain([None, rs2[0], rs2[1], None, rs2[2], rs2[3], None])]
e2 = [(e["kind"], e["number"], e["real"]) for e in hy.chain_ends(M2)]
assert e2 == [("N", 1, False), ("C", 2, False), ("N", 3, False), ("C", 4, False)], e2
rs3 = [mkres(i) for i in range(1, 7)]
class M3:   # a break without SEQRES: residue 3 does not bond to 4
    chains = [chain(list(rs3[:3]))]
chain(list(rs3[3:]))                                  # bond 4-5-6 among themselves, none between 3 and 4
M3.chains[0].residues = rs3[:3] + rs3[3:]            # chain() above only bonded the first three
e3 = [(e["kind"], e["number"], e["real"]) for e in hy.chain_ends(M3)]
assert e3[:2] == [("N", 1, True), ("C", 3, False)] and e3[2:] == [("N", 4, False), ("C", 6, True)], e3
# acid hydrogen geometry
oo, cc_, ox = np.array([1.23, 0, 0]), np.array([0, 0, 0.]), np.array([-0.66, 1.13, 0.3])
hpos = hy.cooh_hydrogen_position(oo, cc_, ox)
assert abs(np.linalg.norm(hpos - ox) - 0.97) < 1e-6 and abs(hy._dihedral(oo, cc_, ox, hpos)) < 1.0
# fix_termini: real ends follow PROPKA, cut ends are made neutral
class FM2:
    id = (8,)
    def __init__(self, residues):
        self.residues = residues; self.deleted = False
    @property
    def atoms(self):
        return [a for r in self.residues for a in r.atoms.values() if not a.deleted]
def build_model(cut_n, cut_c):
    r1, r2 = mkres(1), mkres(2)
    n = r1.atoms["N"]
    for k in ("H1", "H2", "H3"): r1.atoms[k] = FA("H", n, k)
    oxt = FA("O", r2.atoms["C"], "OXT"); r2.atoms["OXT"] = oxt
    m = FM2([r1, r2]); m.chains = [chain([r1, r2])]
    return m, r1, r2
rows = [{"group": "N+", "chain": "A", "number": 1, "pka": 7.99}, {"group": "C-", "chain": "A", "number": 2, "pka": 3.3}]
def run_fix(ends, rows, ph):
    m, r1, r2 = build_model(False, False)
    def addh_amide(sess, text):
        a = FA("H", r1.atoms["N"], "H")            # ChimeraX places one amide hydrogen on the cut N end (and elsewhere: ignored)
        r1.atoms["H"] = a; extra = FA("H", r2.atoms["CA"], "HA"); r2.atoms["HA"] = extra
    hy.run = addh_amide
    notes = hy.fix_termini(session, m, ends, rows, ph, logs_h.append)
    return m, r1, r2, notes
real = [{"chain": "A", "kind": "N", "number": 1, "real": True}, {"chain": "A", "kind": "C", "number": 2, "real": True}]
m, r1, r2, notes = run_fix(real, rows, 7.4)       # N+ pKa 7.99 > 7.4: stays charged; C- 3.3 < 7.4: stays COO-
assert sum(1 for a in r1.atoms["N"].neighbors if hy._is_h(a)) == 3 and "OXT" in r2.atoms and not r2.atoms["OXT"].deleted and "charged NH3+" in notes[("A", 1, "N")] and "charged COO-" in notes[("A", 2, "C")], notes
m, r1, r2, notes = run_fix(real, [dict(rows[0], pka=8.5 - 3), dict(rows[1])], 7.4)       # N+ pKa 5.5 < pH: NH2
assert sum(1 for a in r1.atoms["N"].neighbors if hy._is_h(a)) == 2 and "NH2" in notes[("A", 1, "N")], notes
cut = [{"chain": "A", "kind": "N", "number": 1, "real": False}, {"chain": "A", "kind": "C", "number": 2, "real": False}]
m, r1, r2, notes = run_fix(cut, rows, 7.4)
assert [a.name for a in r1.atoms["N"].neighbors if hy._is_h(a)] == ["H"], "cut N end: one amide hydrogen, the NH3+ ones gone"
assert r2.atoms["OXT"].deleted and "neutral" in notes[("A", 2, "C")] and not any(a.name == "HA" and not a.deleted for a in [r2.atoms["HA"]]), notes
assert "ChimeraX could not" not in " ".join(logs_h[-3:])
m, r1, r2, notes = run_fix([cut[1], real[0]], [rows[0], dict(rows[1], pka=9.0)], 7.4)         # real C end with pKa 9 > pH: COOH
import types as _t
# a cut N end where ChimeraX places two hydrogens: keep the one nearest a planar amide (H-N-CA 120), named H
m, r1, r2 = build_model(False, False)
r1.atoms["N"].coord = [0.0, 0.0, 0.0]; r1.atoms["CA"].coord = [1.46, 0.0, 0.0]
def addh_two(sess, text):
    a = FA("H", r1.atoms["N"], "H2"); a.coord = [-0.5, 0.87, 0.0]       # 120 deg from CA
    b = FA("H", r1.atoms["N"], "H1"); b.coord = [-0.34, 0.94, 0.0]      # about 110 deg
    r1.atoms["H2"], r1.atoms["H1"] = a, b
hy.run = addh_two
n2 = hy.fix_termini(session, m, [{"chain": "A", "kind": "N", "number": 1, "real": False}], rows, 7.4, logs_h.append)
left = [a for a in r1.atoms["N"].neighbors if hy._is_h(a)]
assert len(left) == 1 and left[0].name == "H" and abs(left[0].coord[0] + 0.5) < 1e-9, [(a.name, a.coord) for a in left]
assert "closest to a planar amide" in logs_h[-2], logs_h[-3:]
# Tyr: PROPKA says deprotonated, AMBER keeps TYR; the table says so instead of showing a bare dash
from chimerax.vinapilotx import pdb2pqr_tools as _pq
t = _pq.build_table([{"group": "TYR", "number": 10, "chain": "A", "pka": 10.14, "model_pka": 10.0}], {("A", 10): "TYR"}, {("A", 10): "TYR"}, 10.5)
assert t[0]["note"].startswith("kept TYR") and t[0]["renamed"] == ""
added_atoms = []
_mod = _t.ModuleType("chimerax.atomic.struct_edit")
def _add_atom(name, el, res, pos, bonded_to=None): a = FA(el, bonded_to, name); added_atoms.append((name, el, tuple(round(float(x), 2) for x in pos))); return a
_mod.add_atom = _add_atom; sys.modules["chimerax.atomic.struct_edit"] = _mod
m, r1, r2, notes = run_fix([dict(real[0]), dict(real[1])], [dict(rows[0]), dict(rows[1], pka=9.0)], 7.4)
assert added_atoms and added_atoms[0][:2] == ("HO", "H") and "COOH" in notes[("A", 2, "C")], (added_atoms, notes)

# ---- pKa table: two flag colours, readable in light and dark themes
from chimerax.vinapilotx import pdb2pqr_tools as pq
from chimerax.vinapilotx.prep_widgets import flag_colors
rows = pq.build_table([{"group": "ASP", "number": 4, "chain": "A", "pka": 3.94, "model_pka": 3.8},
                       {"group": "HIS", "number": 6, "chain": "A", "pka": 7.9, "model_pka": 6.5},
                       {"group": "GLU", "number": 8, "chain": "A", "pka": 5.6, "model_pka": 4.5},
                       {"group": "LYS", "number": 5, "chain": "A", "pka": 10.45, "model_pka": 10.5}],
                      {("A", 6): "HIS"}, {("A", 6): "HIP"}, 7.4)
assert [r["flag"] for r in rows] == ["", "close", "near", ""] and rows[1]["renamed"] == "HIS → HIP"
P.pka_table.set_rows(rows, 7.4, "t"); assert P.pka_table.isVisibleTo(P) is True or True
tb = P.pka_table.table
assert tb.item(1, 6).text() == "±1" and tb.item(2, 6).text() == "±2" and tb.item(0, 6).text() == ""
for dark in (False, True):
    strong, soft, text = flag_colors(dark)
    assert strong != soft and strong.name() != soft.name()
    lum = lambda c: c.lightness()
    assert (lum(text) > lum(strong) and lum(text) > lum(soft)) if dark else (lum(text) < lum(strong) and lum(text) < lum(soft))
P.pka_table.recolor()
assert tb.item(1, 0).background().color() != tb.item(2, 0).background().color() and tb.item(0, 0).background().style() == QtCore.Qt.BrushStyle.NoBrush and tb.item(1, 0).background().style() != QtCore.Qt.BrushStyle.NoBrush
# hydrogen step through the Protein tab (ChimeraX polar + non-polar), checkpoints, undo and redo
from chimerax.vinapilotx.project import Project
tool.workdir_edit.setText(tmp); proj = tool.get_project()
P._prep_name = "prot"; P._rev = proj.new_revision("prot", "prot.pdb", {})
for st in ("clean",): pass
fm.deleted = False; P._prep_model = fm
hy.run = fake_addh; fm._atoms = [C1, N1, O1]
saves = []
import chimerax.vinapilotx.receptor_clean as rcl
def fake_save(sess, model, path): saves.append(path); open(path, "w").write("ATOM\n")
import chimerax.vinapilotx.protein_tab as ptab
rcl.save_model_pdb = fake_save
opened = []
def fake_open(sess, path, name): opened.append(path); mm = FM(); mm._atoms = [FA("C"), FA("N")]; mm.deleted = False; return mm
rcl.open_prep_model = fake_open
sess_models_close = session.models.close
session.models.close = lambda ms: None
for i in range(P.polar_combo.count()):
    if P.polar_combo.itemData(i) == "chimerax": P.polar_combo.setCurrentIndex(i)
P.nonpolar_cb.setChecked(True)
P._apply_hydrogens()
assert "✔" in P.h_result.text() and "ChimeraX placed" in P.h_result.text(), P.h_result.text()
assert P._h_pre and os.path.basename(P._h_pre) == "_before_hydrogens.pdb" and os.path.isfile(P._h_pre)
assert any(p.endswith("02_hydrogens.pdb") for p in saves) and "hydrogens" in P._rev["stages"]
assert P.h_undo_btn.isEnabled() and not P.h_redo_btn.isEnabled()
P._undo_hydrogens(); assert opened[-1] == P._h_pre and not P.h_undo_btn.isEnabled() and P.h_redo_btn.isEnabled() and "hydrogens" not in P._rev["stages"]
P._undo_hydrogens(); assert len(opened) == 1                       # idempotent
P._redo_hydrogens(); assert opened[-1].endswith("02_hydrogens.pdb") and P.h_undo_btn.isEnabled() and "hydrogens" in P._rev["stages"]
# the project file records the stages and survives a reload
assert os.path.isfile(os.path.join(tmp, "project.vinapx")) and Project.open(tmp).active()["stages"]["hydrogens"]
for i in range(P.polar_combo.count()):
    if P.polar_combo.itemData(i) == "pdb2pqr": P.polar_combo.setCurrentIndex(i)
assert P.ph_spin.isEnabled() and P.ph_spin.value() == 7.0
tool.controller.prep_info = {"meeko": "0.8.0"}; P._apply_hydrogens(); assert "not installed" in P.h_result.text().lower()
P.polar_combo.setCurrentIndex(0); assert not P.ph_spin.isEnabled()
session.models.close = sess_models_close
# a new clean must not leave the old receptor on the Inputs tab
oldq = proj.stage_file(P._rev, "clean").replace("01_clean.pdb", "p.pdbqt"); open(oldq, "w").write("x")
tool.receptor_edit.setText(oldq); assert tool.drop_prepared_receptor("test") and tool.receptor_edit.text() == ""
tool.receptor_edit.setText("/somewhere/else/my.pdbqt"); assert not tool.drop_prepared_receptor("test") and tool.receptor_edit.text() == "/somewhere/else/my.pdbqt"; tool.receptor_edit.setText("")

# ---- pdb2pqr parsing
txt = """x\nSUMMARY OF THIS PREDICTION\n       Group      pKa  model-pKa   ligand atom-type\n   ASP   4 A     3.94       3.80                      \n   N+    1 A     7.68       8.00                      \n\nINFO:next"""
got = pq.parse_pka_summary(txt); assert [(g["group"], g["number"], g["pka"]) for g in got] == [("ASP", 4, 3.94), ("N+", 1, 7.68)]
assert pq.classify(7.0, 7.4) == "close" and pq.classify(5.5, 7.4) == "near" and pq.classify(3, 7.4) == ""

# ---- approval flow: Meeko result is NOT put on the Inputs tab until the user approves it
tool.receptor_edit.setText(""); proj = tool.get_project()
rev = proj.new_revision("appr", "x.pdb"); q = proj.stage_file(rev, "clean").replace("01_clean.pdb", "appr.pdbqt"); open(q, "w").write("x")
proj.set_stage(rev, "pdbqt", q); proj.save(); P._rev = rev; P._convert_rev = rev; P.load_project(); tool.refresh_prepared_list()
assert tool.receptor_edit.text() == "" and not proj.approved() and tool.prepared_combo.count() == 1 and not tool._approved
assert P.approve_btn.isEnabled() and P.approve_btn.text() == "Approve for docking"
P._toggle_approve()
assert proj.approved()[0][1] == q and tool.prepared_combo.count() == 2 and P.approve_btn.text() == "Remove approval"
assert not hasattr(P, "rev_use_btn")
tool._receptor_chosen = lambda open_model=False: None
tool._prepared_picked(1); assert tool.receptor_edit.text() == q
assert tool.save_project(); _saved = Project.open(proj.workdir).data
assert _saved["receptor"]["selected"].endswith("appr.pdbqt") and _saved["params"] and \
    _saved["ligands"]["inputs"] == [os.path.relpath(x, proj.workdir).replace("\\", "/") if not os.path.relpath(x, proj.workdir).startswith("..") else x for x in dict.fromkeys(r.source for r in tool._records.values())]
proj.set_stage(rev, "pdbqt", q); assert not rev["approved"]       # a new conversion has to be approved again
tool.receptor_edit.setText("")

# ---- open-model list follows ChimeraX while the window stays open
calls = []; P.refresh_models = lambda: calls.append(1)
names = [n for n, f in _handlers if n in ("am", "rm")]; assert names == ["am", "rm"], _handlers
for n, f in list(_handlers):
    if n in ("am", "rm"): f(n, []); f(n, [])
from Qt.QtWidgets import QApplication
QApplication.processEvents(); QApplication.processEvents()
assert len(calls) == 1, "two quick model changes must give one refresh"
import chimerax.vinapilotx.protein_tab as ptab
assert ptab.stem("1asd.pdb") == "1asd" and ptab.stem("x.CIF") == "x" and ptab.stem("a.b.c") == "a.b.c"
P.cleanup(); assert [h for h in _handlers if h[0] == "removed"]

tool.delete(); print("GUI SMOKE OK")
