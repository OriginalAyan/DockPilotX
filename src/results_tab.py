"""Results tab: expandable ligand -> pose tree, ratings, per-pose H-bond/clash analysis, exports."""

import json
import os
import shutil

from Qt.QtCore import Qt, QTimer
from Qt.QtGui import QGuiApplication
from Qt.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QPushButton,
    QSpinBox, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from . import docking, layout
from .analysis_store import AnalysisStore
from .pose_viewer import PoseViewer

COLS = ["Ligand / pose", "Affinity (kcal/mol)", "RMSD l.b.", "RMSD u.b.", "Rating", "H-bonds", "Clashes",
        "Poses", "Seed", "Status"]
C_NAME, C_AFF, C_LB, C_UB, C_RATE, C_HB, C_CL, C_POSES, C_SEED, C_STATUS = range(10)
STARS = ["—", "★", "★★", "★★★", "★★★★", "★★★★★"]


def hb_text(rec):
    if rec.get("hbonds") is None:
        return ""
    n = rec.get("hbonds_relaxed") or 0
    return "%d (+%d)" % (rec["hbonds"], n) if n else str(rec["hbonds"])


class _Item(QTreeWidgetItem):
    """Tree item sorted numerically; poses always stay in pose order inside their ligand."""

    def __init__(self, texts, job=None, pose=None):
        super().__init__(texts)
        self.job = job
        self.pose = pose            # None for a ligand row
        self.sort_values = {}

    def __lt__(self, other):
        tree = self.treeWidget()
        col = tree.sortColumn() if tree is not None else 0
        if self.pose is not None and getattr(other, "pose", None) is not None:
            return self.pose < other.pose
        a, b = self.sort_values.get(col), getattr(other, "sort_values", {}).get(col)
        if a is not None and b is not None:
            return a < b
        if a is None and b is not None:
            return False
        if a is not None and b is None:
            return True
        return self.text(col).lower() < other.text(col).lower()


class ResultsTab(QWidget):
    def __init__(self, session, controller, host):
        """*host*: the tool, providing workdir_edit and the receptor-model helpers used by PoseViewer."""
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host
        self.jobs = []
        self.store = None
        self.workdir = None
        self.viewer = PoseViewer(session, host, max_open=int(controller.settings.get("max_open_files", 5)),
                                colors=controller.settings.get("colors"))
        self._ligand_items = {}
        self._building = False
        self._sel_timer = QTimer(self)
        self._sel_timer.setSingleShot(True)
        self._sel_timer.setInterval(150)
        self._sel_timer.timeout.connect(self._apply_selection)
        self._build_ui()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addWidget(QLabel("Max ligand files kept open"))
        self.max_open_spin = QSpinBox()
        self.max_open_spin.setRange(1, 50)
        self.max_open_spin.setValue(self.viewer.max_open)
        self.max_open_spin.setToolTip("Each ligand file holds all its poses. The oldest file is closed beyond this number.")
        self.max_open_spin.valueChanged.connect(self._max_open_changed)
        top.addWidget(self.max_open_spin)
        self.receptor_cb = QCheckBox("Show receptor")
        self.receptor_cb.setChecked(True)
        top.addWidget(self.receptor_cb)
        top.addStretch(1)
        lay.addLayout(top)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(len(COLS))
        self.tree.setHeaderLabels(COLS)
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.setSortingEnabled(True)
        self.tree.setRootIsDecorated(True)
        self.tree.setUniformRowHeights(False)
        self.tree.setMinimumHeight(240)
        self.tree.header().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.tree.header().setStretchLastSection(True)
        self.tree.setColumnWidth(C_NAME, 160)
        self.tree.itemExpanded.connect(self._expanded)
        self.tree.itemSelectionChanged.connect(lambda: self._sel_timer.start())
        lay.addWidget(self.tree, 1)

        row = QHBoxLayout()
        self.hb_btn = QPushButton("H-bonds")
        self.hb_btn.setToolTip("Calculate H-bonds for the selected pose(s).")
        self.hb_btn.clicked.connect(lambda: self.analyze("hbonds"))
        self.cl_btn = QPushButton("Clashes")
        self.cl_btn.setToolTip("Calculate clashes for the selected pose(s).")
        self.cl_btn.clicked.connect(lambda: self.analyze("clashes"))
        self.vd_btn = QPushButton("Open in ViewDock")
        self.vd_btn.clicked.connect(self.open_viewdock)
        for b in (self.hb_btn, self.cl_btn, self.vd_btn):
            row.addWidget(b)
        lay.addLayout(row)

        row = QHBoxLayout()
        export = QPushButton("Export CSV…")
        export.clicked.connect(self.export_csv)
        row.addWidget(export)
        save = QPushButton("Save rated poses…")
        save.setToolTip("Copies poses at or above this rating (and a CSV) into vinax_selected_poses/.")
        save.clicked.connect(self.save_rated)
        row.addWidget(save)
        self.min_rating = QSpinBox()
        self.min_rating.setRange(1, 5)
        self.min_rating.setValue(4)
        self.min_rating.setPrefix("≥ ")
        self.min_rating.setSuffix(" ★")
        row.addWidget(self.min_rating)
        load = QPushButton("Load folder…")
        load.clicked.connect(self.load_folder)
        row.addWidget(load)
        lay.addLayout(row)

        self.note = QLabel("Select a pose to show it. Expand a ligand to see all poses.")
        self.note.setWordWrap(True)
        lay.addWidget(self.note)

    def _max_open_changed(self, v):
        self.viewer.max_open = v
        self.controller.settings.set("max_open_files", v)

    # ------------------------------------------------------------ populate
    def set_jobs(self, jobs, workdir, receptor_file):
        """Show *jobs* (may still be running) for *workdir*."""
        self.viewer.set_project(workdir, receptor_file)
        self.workdir = workdir
        self.store = AnalysisStore(workdir)
        self.jobs = list(jobs)
        self._building = True
        self.tree.setSortingEnabled(False)
        self.tree.clear()
        self._ligand_items = {}
        for job in self.jobs:
            item = _Item([job.name], job=job)
            item.addChild(QTreeWidgetItem(["…"]))  # placeholder so the expand arrow appears
            self.tree.addTopLevelItem(item)
            self._ligand_items[job.name] = item
            self._fill_ligand_row(item)
        self.tree.setSortingEnabled(True)
        self._building = False

    def _best_rating(self, job):
        if self.store is None:
            return 0
        return max([int(self.store.get(job.name, n).get("rating", 0) or 0) for n in range(1, job.n_poses + 1)],
                   default=0)

    def _fill_ligand_row(self, item):
        job = item.job
        item.setText(C_AFF, "" if job.best is None else "%.3f" % job.best)
        item.sort_values[C_AFF] = job.best
        item.setText(C_POSES, str(job.n_poses) if job.n_poses else "")
        item.sort_values[C_POSES] = job.n_poses
        item.setText(C_SEED, "" if job.seed is None else str(job.seed))
        item.sort_values[C_SEED] = job.seed
        r = self._best_rating(job)
        item.setText(C_RATE, STARS[r] if r else "")
        item.sort_values[C_RATE] = r
        status = job.status
        if job.status == docking.STATUS_RUNNING:
            status = "running %d%%" % job.progress
        elif job.error:
            status = "%s: %s" % (job.status, job.error)
        item.setText(C_STATUS, status)
        item.setToolTip(C_STATUS, job.error or "")

    def update_job(self, job):
        item = self._ligand_items.get(job.name)
        if item is None:
            return
        item.job = job
        self._fill_ligand_row(item)
        if item.isExpanded():
            self._fill_children(item)

    def _expanded(self, item):
        if getattr(item, "pose", None) is None and getattr(item, "job", None) is not None:
            self._fill_children(item)

    def _fill_children(self, item):
        job = item.job
        self._building = True
        item.takeChildren()
        for n, (aff, lb, ub) in enumerate(job.poses, start=1):
            rec = self.store.get(job.name, n) if self.store else {}
            child = _Item(["Pose %d" % n, "%.3f" % aff, "%.3f" % lb, "%.3f" % ub, "",
                           hb_text(rec),
                           "" if rec.get("clashes") is None else str(rec["clashes"])], job=job, pose=n)
            for col, v in ((C_AFF, aff), (C_LB, lb), (C_UB, ub), (C_HB, (rec.get("hbonds") or 0) + (rec.get("hbonds_relaxed") or 0) if rec.get("hbonds") is not None else None), (C_CL, rec.get("clashes"))):
                child.sort_values[col] = v
            tips = rec.get("hbond_partners")
            if tips:
                child.setToolTip(C_HB, "\n".join(tips))
            item.addChild(child)
            combo = QComboBox()
            combo.addItems(STARS)
            combo.setCurrentIndex(int(rec.get("rating", 0) or 0))
            combo.currentIndexChanged.connect(lambda idx, j=job, p=n: self._rate(j, p, idx))
            self.tree.setItemWidget(child, C_RATE, combo)
        self._building = False

    # ----------------------------------------------------------- selection
    def _selected_poses(self):
        out, seen = [], set()
        for it in self.tree.selectedItems():
            job = getattr(it, "job", None)
            if job is None or not job.poses:
                continue
            pose = it.pose if it.pose is not None else 1  # a ligand row stands for its best pose
            if (job.name, pose) not in seen:
                seen.add((job.name, pose))
                out.append((job, pose))
        return out

    def _apply_selection(self):
        if self._building or self.workdir is None:
            return
        wanted = self._selected_poses()
        try:
            if wanted:
                if len({j.name for j, _ in wanted}) > self.viewer.max_open:
                    self.note.setText("Showing the first %d ligands (see 'Max ligand files kept open')."
                                      % self.viewer.max_open)
                else:
                    self.note.setText("Select a pose to show it. Expand a ligand to see all poses.")
                self.viewer.show(wanted, with_receptor=self.receptor_cb.isChecked())
            else:
                self.viewer.hide_all()
        except Exception as err:
            self.note.setText("Could not show the pose: %s" % err)
            self.session.logger.warning("Could not show pose: %s" % err)

    # ------------------------------------------------------------- ratings
    def _rate(self, job, pose, idx):
        if self._building or self.store is None:
            return
        self.store.set(job.name, pose, rating=idx)
        item = self._ligand_items.get(job.name)
        if item is not None:
            self._fill_ligand_row(item)

    # ------------------------------------------------------------ analysis
    def analyze(self, kind):
        poses = [(it.job, it.pose) for it in self.tree.selectedItems() if it.pose is not None]
        if not poses:
            self.note.setText("Select one or more poses (expand a ligand first).")
            return
        if len(poses) > 3:
            self.note.setText("Calculating for the first 3 of %d selected poses; repeat for the others." % len(poses))
            poses = poses[:3]
        label = "H-bonds" if kind == "hbonds" else "clashes"
        QGuiApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            for job, pose in poses:
                self.note.setText("Calculating %s for %s pose %d …" % (label, job.name, pose))
                QGuiApplication.processEvents()
                try:
                    count, relaxed, partners = self.viewer.analyze(job, pose, kind)
                except Exception as err:
                    self.note.setText("%s failed: %s" % (label, err))
                    self.session.logger.warning("%s calculation failed: %s" % (label, err))
                    return
                if kind == "hbonds":
                    self.store.set(job.name, pose, hbonds=count, hbonds_relaxed=relaxed, hbond_partners=partners)
                else:
                    self.store.set(job.name, pose, clashes=count)
                self._refresh_pose_row(job, pose)
            self.note.setText("Done. The %s are saved and shown whenever the pose is shown." % label)
        finally:
            QGuiApplication.restoreOverrideCursor()

    def _refresh_pose_row(self, job, pose):
        item = self._ligand_items.get(job.name)
        if item is None:
            return
        for i in range(item.childCount()):
            ch = item.child(i)
            if ch.pose == pose:
                rec = self.store.get(job.name, pose)
                ch.setText(C_HB, hb_text(rec))
                ch.setText(C_CL, "" if rec.get("clashes") is None else str(rec["clashes"]))
                ch.sort_values[C_HB] = None if rec.get("hbonds") is None else rec["hbonds"] + (rec.get("hbonds_relaxed") or 0)
                ch.sort_values[C_CL] = rec.get("clashes")
                if rec.get("hbond_partners"):
                    ch.setToolTip(C_HB, "\n".join(rec["hbond_partners"]))

    # ------------------------------------------------------------- viewdock
    def open_viewdock(self):
        from chimerax.core.commands import run
        done = set()
        for it in self.tree.selectedItems():
            job = getattr(it, "job", None)
            if job is None or not job.poses or job.name in done or not os.path.isfile(job.out_path):
                continue
            done.add(job.name)
            try:
                run(self.session, "viewdock %s" % self.viewer.ligand_spec(job))
            except Exception as err:
                self.note.setText("ViewDock failed: %s" % err)

    # --------------------------------------------------------------- export
    def export_csv(self):
        if not self.jobs:
            return
        start = os.path.join(self.workdir or "", "results_all_poses.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export results", start, "CSV (*.csv)")
        if not path:
            return
        docking.write_pose_csv(self.jobs, path, self.store)
        self.note.setText("Exported all poses to %s" % path)

    def save_rated(self):
        if not self.jobs or self.store is None:
            return
        minimum = self.min_rating.value()
        by_name = {j.name: j for j in self.jobs}
        chosen = [(n, p) for n, p, _ in self.store.rated(minimum) if n in by_name and by_name[n].poses]
        if not chosen:
            self.note.setText("No poses are rated %d★ or higher." % minimum)
            return
        out_dir = layout.make(self.workdir, "selected")
        for n, p in chosen:
            docking.extract_pose(by_name[n].out_path, p, os.path.join(out_dir, "%s_pose%d.pdbqt" % (n, p)))
        wanted = set(chosen)
        docking.write_pose_csv(self.jobs, os.path.join(out_dir, "selected_poses.csv"), self.store,
                               only=lambda name, pose, rec: (name, pose) in wanted)
        rec_file = self.viewer.receptor_file
        if rec_file and os.path.isfile(rec_file):
            shutil.copyfile(rec_file, os.path.join(out_dir, os.path.basename(rec_file)))
        self.note.setText("Saved %d pose(s) rated ≥ %d★ to %s" % (len(chosen), minimum, out_dir))

    # ----------------------------------------------------------------- load
    def load_folder(self):
        start = self.host.workdir_edit.text().strip() if hasattr(self.host, "workdir_edit") else ""
        d = QFileDialog.getExistingDirectory(self, "Select a working directory", start)
        if d:
            self.load_workdir(d)

    def load_workdir(self, d):
        jobs = docking.load_results_dir(d)
        if not jobs:
            self.note.setText("No Vina results found in %s" % layout.existing(d, "results"))
            return False
        receptor = None
        try:
            with open(os.path.join(d, "manifest.json"), encoding="utf-8") as f:
                used = json.load(f).get("receptor_used")
            if used and os.path.isfile(used):
                receptor = used
        except (OSError, ValueError):
            pass
        if receptor is None:
            rdir = layout.existing(d, "receptor")
            if os.path.isdir(rdir):
                for f in sorted(os.listdir(rdir)):
                    if f.lower().endswith(".pdbqt"):
                        receptor = os.path.join(rdir, f)
                        break
        if hasattr(self.host, "workdir_edit"):
            self.host.workdir_edit.setText(d)
        self.set_jobs(jobs, d, receptor)
        n_rated = len(self.store.rated(1))
        self.note.setText("Loaded %d ligand(s) from %s%s" % (
            len(jobs), d, " (restored %d rated pose(s) and saved analysis)" % n_rated if n_rated else ""))
        return True

    def sort_by_best(self):
        self.tree.sortItems(C_AFF, Qt.SortOrder.AscendingOrder)

    def close(self):
        try:
            self.viewer.close_all()
        except Exception:
            pass
