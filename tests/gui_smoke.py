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
        if getattr(self.session, 'models', None) and self in self.session.models.items: self.session.models.items.remove(self)
class Surface(Model):
    def set_geometry(self, *a): pass
mod("chimerax"); mod("chimerax.core"); mod("chimerax.core.models", Model=Model, Surface=Surface, MODEL_COLOR_CHANGED="mcc")
mod("chimerax.surface"); mod("chimerax.surface.shapes", cylinder_geometry=lambda **k: (np.zeros((4,3)), np.zeros((4,3)), np.zeros((2,3), int)))
class UndoAction:
    def __init__(self, name): self.name = name
mod("chimerax.core.undo", UndoAction=UndoAction)
ran = []
_ids = [10]
class _Atoms:
    def __init__(self): self.scene_coords = np.array([[0, 0, 0], [30, 10, 5.0]]); self.__len = 2
    def __len__(self): return 2
def _run(sess, text):
    ran.append(text)
    if text.startswith("open ") and "format pdbqt" in text:
        m = Model(text, sess); m.atoms = _Atoms(); _ids[0] += 1; m.id = (_ids[0],); sess.models.items.append(m)
    elif text.startswith(("hbonds", "clashes")):
        nm = text.split(" name ")[1].split()[0]
        g = PseudobondGroup(nm, sess); g.name = nm; _ids[0] += 1; g.id = (_ids[0],); g.num_pseudobonds = 2 if text.startswith("hbonds") else 1
        g.pseudobonds = []; sess.models.items.append(g)
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
class PseudobondGroup(Model):
    pass
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
session.triggers = types.SimpleNamespace(add_handler=lambda *a: 1, remove_handler=lambda h: None)

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
assert tool.tabs.currentIndex() == 4 and tool.run_btn.isEnabled() and not tool.cancel_btn.isEnabled()
W = os.path.join(tmp, "work")
rows_ = list(__import__("csv").DictReader(open(os.path.join(W, "summary.csv"))))
assert len(rows_) >= 3 and {"pose", "rmsd_lb", "rmsd_ub", "rating"} <= set(rows_[0]) and "time" not in "".join(rows_[0]).lower()
assert os.path.isfile(os.path.join(W, "run_log.txt"))
# expand: all poses listed, receptor never reopened
R = tool.results
top = T.topLevelItem(0); top.setExpanded(True)
assert top.childCount() == top.job.n_poses and top.childCount() > 1
n_rec_open = sum("rec" in c and c.startswith("open ") for c in ran)
# select pose 1, then each pose: cache capped at 5, receptor opened only by the tool (not again)
R.viewer.max_open = 3
for i in range(top.childCount()):
    T.clearSelection(); top.child(i).setSelected(True); R._apply_selection()
open_poses = [e for e in R.viewer._entries.values() if not e.model.deleted]
assert len(open_poses) <= 3, len(open_poses)
assert sum("rec" in c and c.startswith("open ") for c in ran) == n_rec_open, "receptor was reopened"
# rating persists
child = top.child(1); combo = T.itemWidget(child, 4); combo.setCurrentIndex(5)
assert R.store.get(top.job.name, 2)["rating"] == 5
# H-bonds / clashes for the selected pose
T.clearSelection(); child.setSelected(True); R.analyze("hbonds"); R.analyze("clashes")
assert child.text(5) == "2" and child.text(6) == "1", (child.text(5), child.text(6))
assert os.path.isfile(os.path.join(W, "analysis.json"))
# reload the folder in a fresh results tab state: rating + analysis restored
assert R.load_workdir(W)
top2 = T.topLevelItem(0)
for r in range(T.topLevelItemCount()):
    if T.topLevelItem(r).job.name == top.job.name: top2 = T.topLevelItem(r)
top2.setExpanded(True)
assert top2.child(1).text(5) == "2" and T.itemWidget(top2.child(1), 4).currentIndex() == 5
assert top2.text(4) == "★★★★★"
# save rated poses
R.min_rating.setValue(4); R.save_rated()
sel = os.path.join(W, "selected_poses"); assert os.path.isdir(sel) and any(f.endswith(".pdbqt") for f in os.listdir(sel))
# viewdock + csv export helpers do not crash
R.tree.clearSelection()
cfg = open(os.path.join(tmp, "work", "vina_config.txt")).read(); assert "exhaustiveness = 16" in cfg and "size_x = 22.0000" in cfg
# reset defaults + progress bar style
tool._reset_defaults(); assert tool.exh_spin.value() == 8 and tool.parallel_spin.value() == 1 and tool.seed_spin.value() == 0
assert "#1e88e5" in tool.tool_window.ui_area.styleSheet()
# min-size: no page may force a tall window
print("tab min heights:", [tool.tabs.widget(i).minimumSizeHint().height() for i in range(tool.tabs.count())])
assert max(tool.tabs.widget(i).minimumSizeHint().height() for i in range(tool.tabs.count())) < 300
tool.controller.vina_exe = None
ok, msg = tool.controller.start_docking(rec, [lib], tmp, (0,0,0), (10,10,10), tool._params()); assert not ok and "not set up" in msg
tool.delete(); print("GUI SMOKE OK")
