"""The Lattice VinaPilotX tool window: Settings, Inputs, Box, Run and Results tabs."""

import csv
import os

from Qt.QtCore import Qt, QUrl
from Qt.QtGui import QBrush, QColor, QDesktopServices, QFont, QGuiApplication
from Qt.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QFrame,
    QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QPlainTextEdit,
    QProgressBar, QPushButton, QScrollArea, QSpinBox, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)
from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run

from . import docking, pdbqt_tools
from .box_panel import BoxPanel
from .controller import get_controller

OK_STYLE = "color: #2e7d32; font-weight: bold;"
BAD_STYLE = "color: #c62828; font-weight: bold;"
WARN_STYLE = "color: #b26a00;"
HIGHLIGHT = "#1f4e6b"  # dark blue: readable with white text on ChimeraX's dark theme

# Blue progress bars on a light track, so empty and filled are clearly different.
TOOL_STYLE = """
QProgressBar { border: 1px solid #8a8a8a; border-radius: 3px; background: #e9eef3;
               text-align: center; color: #111111; min-height: 16px; }
QProgressBar::chunk { background-color: #1e88e5; border-radius: 2px; }
"""

DEFAULTS = {"exhaustiveness": 8, "num_modes": 9, "energy_range": 3.0, "scoring": "vina",
            "cpu": 0, "parallel_jobs": 1, "seed": 0}


class _NumItem(QTableWidgetItem):
    """Table item that sorts by a numeric value (missing values last)."""

    def __init__(self, text, value=None):
        super().__init__(text)
        self._value = value

    def __lt__(self, other):
        a = self._value if self._value is not None else float("inf")
        b = getattr(other, "_value", None)
        b = b if b is not None else float("inf")
        return a < b


def _row_with_browse(line_edit, button_text, slot):
    row = QHBoxLayout()
    row.addWidget(line_edit, 1)
    b = QPushButton(button_text)
    b.clicked.connect(slot)
    row.addWidget(b)
    return row, b


def _scrolled(widget):
    """Wrap a tab page in a scroll area so the tool never forces a tall minimum window height."""
    area = QScrollArea()
    area.setWidgetResizable(True)
    area.setFrameShape(QFrame.Shape.NoFrame)
    area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
    area.setWidget(widget)
    return area


def _chimerax_path(path):
    return os.path.abspath(path).replace("\\", "/")


def _fmt_dims(dims):
    return " × ".join("%.1f" % d for d in dims)


class VinaPilotXTool(ToolInstance):
    SESSION_ENDURING = False
    help = "help:user/tools/vinapilotx.html"

    def __init__(self, session, tool_name):
        super().__init__(session, tool_name)
        self.display_name = "Lattice VinaPilotX"
        from chimerax.ui import MainToolWindow
        self.tool_window = MainToolWindow(self)
        self.controller = get_controller(session)
        self._rows = {}               # results: job name -> (items, job)
        self._records = {}            # ligand table: id -> LigandRecord
        self._used_names = set()
        self._next_rid = 0
        self._slots = []              # per-running-ligand progress widgets
        self._receptor_path = None
        self._receptor_models = []
        self._build_ui()
        self._restore_settings()
        self.controller.add_listener(self._on_event)
        self._refresh_vina_status()
        self.tool_window.manage(placement="side")

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        parent = self.tool_window.ui_area
        parent.setStyleSheet(TOOL_STYLE)
        parent.setMinimumSize(300, 200)  # small on purpose: pages scroll instead of forcing a tall window
        outer = QVBoxLayout(parent)
        outer.setContentsMargins(4, 4, 4, 4)
        self.tabs = QTabWidget()
        outer.addWidget(self.tabs)
        self.tabs.addTab(_scrolled(self._build_settings_tab()), "1 Settings")
        self.tabs.addTab(_scrolled(self._build_inputs_tab()), "2 Inputs")
        self.box_panel = BoxPanel(self.session)
        self.tabs.addTab(self.box_panel, "3 Box")
        self.tabs.addTab(_scrolled(self._build_run_tab()), "4 Run")
        self.tabs.addTab(_scrolled(self._build_results_tab()), "5 Results")
        self.tabs.currentChanged.connect(self._tab_changed)

    # ---- Settings ----------------------------------------------------
    def _build_settings_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox("AutoDock Vina program")
        gl = QVBoxLayout(g)
        self.vina_status = QLabel("")
        self.vina_status.setWordWrap(True)
        gl.addWidget(self.vina_status)
        self.vina_edit = QLineEdit()
        self.vina_edit.setPlaceholderText("Path to vina / vina.exe, or the folder that contains it")
        self.vina_edit.editingFinished.connect(self._verify_vina)
        row, _ = _row_with_browse(self.vina_edit, "Browse…", self._browse_vina)
        gl.addLayout(row)
        btns = QHBoxLayout()
        self.verify_btn = QPushButton("Verify")
        self.verify_btn.clicked.connect(self._verify_vina)
        self.download_btn = QPushButton("Download latest Vina")
        self.download_btn.setToolTip("Download the official Vina binary from GitHub releases.")
        self.download_btn.clicked.connect(self.controller.download_vina)
        self.cancel_dl_btn = QPushButton("Cancel")
        self.cancel_dl_btn.clicked.connect(self.controller.cancel_download)
        self.cancel_dl_btn.setVisible(False)
        for b in (self.verify_btn, self.download_btn, self.cancel_dl_btn):
            btns.addWidget(b)
        btns.addStretch(1)
        gl.addLayout(btns)
        self.dl_bar = QProgressBar()
        self.dl_bar.setVisible(False)
        gl.addWidget(self.dl_bar)
        note = QLabel("Vina runs as a separate program. Nothing is installed into ChimeraX.")
        note.setWordWrap(True)
        gl.addWidget(note)
        lay.addWidget(g)

        g = QGroupBox("Working directory")
        gl = QVBoxLayout(g)
        self.workdir_edit = QLineEdit()
        self.workdir_edit.setPlaceholderText("Choose or type a folder, e.g. C:/docking/project1")
        row, _ = _row_with_browse(self.workdir_edit, "Browse…", self._browse_workdir)
        gl.addLayout(row)
        lay.addWidget(g)

        g = QGroupBox("Opening .pdbqt files")
        gl = QVBoxLayout(g)
        self.native_cb = QCheckBox("Open PDBQT files automatically")
        self.native_cb.setToolTip("Forces ChimeraX's built-in reader for .pdbqt files and ViewDock multi-pose "
                                  "rendering. Uncheck if file-opening conflicts occur with other extensions.")
        self.native_cb.toggled.connect(self._native_toggled)
        gl.addWidget(self.native_cb)
        n = QLabel("Applies to file opening and drag-and-drop. Uncheck to restore default behavior.")
        n.setWordWrap(True)
        gl.addWidget(n)
        lay.addWidget(g)
        lay.addStretch(1)
        return w

    def _native_toggled(self, checked):
        if getattr(self, "_restoring", False):
            return
        ok = self.controller.set_native_pdbqt(checked)
        if not ok:
            self.native_cb.blockSignals(True)
            self.native_cb.setChecked(not checked)
            self.native_cb.blockSignals(False)

    def _browse_vina(self):
        path, _ = QFileDialog.getOpenFileName(self.tool_window.ui_area, "Select the AutoDock Vina program",
                                              self.vina_edit.text() or "")
        if path:
            self.vina_edit.setText(path)
            self._verify_vina()

    def _verify_vina(self):
        text = self.vina_edit.text().strip()
        if text:
            self.controller.set_vina(text)

    def _refresh_vina_status(self):
        c = self.controller
        if c.vina_exe:
            self.vina_status.setText("✔ AutoDock Vina %s is ready." % c.vina_version)
            self.vina_status.setStyleSheet(OK_STYLE)
            self.vina_edit.setText(c.vina_exe)
        else:
            self.vina_status.setText("✘ AutoDock Vina is not set up yet. Browse to vina.exe, or press "
                                     "“Download latest Vina”.")
            self.vina_status.setStyleSheet(BAD_STYLE)

    def _browse_workdir(self):
        d = QFileDialog.getExistingDirectory(self.tool_window.ui_area, "Select the working directory",
                                             self._start_dir(self.workdir_edit.text()))
        if d:
            self.workdir_edit.setText(d)

    # ---- Inputs ------------------------------------------------------
    def _build_inputs_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)

        g = QGroupBox("Receptor (prepared .pdbqt)")
        gl = QVBoxLayout(g)
        self.receptor_edit = QLineEdit()
        self.receptor_edit.setPlaceholderText("Prepared receptor PDBQT file")
        self.receptor_edit.editingFinished.connect(lambda: self._receptor_chosen(open_model=True))
        row, _ = _row_with_browse(self.receptor_edit, "Browse…", self._browse_receptor)
        gl.addLayout(row)
        self.receptor_status = QLabel("")
        self.receptor_status.setWordWrap(True)
        gl.addWidget(self.receptor_status)
        lay.addWidget(g)

        g = QGroupBox("Ligands (.pdbqt)")
        gl = QVBoxLayout(g)
        row = QHBoxLayout()
        for text, slot in (("Add files…", self._add_ligand_files), ("Add folder…", self._add_ligand_folder),
                           ("Remove selected", self._remove_ligands), ("Clear", self._clear_ligands),
                           ("Save table (CSV)", self._save_ligand_csv)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        gl.addLayout(row)
        self.recursive_cb = QCheckBox("Include sub-folders")
        gl.addWidget(self.recursive_cb)
        cols = ["Ligand", "Heavy atoms", "Torsions", "Dimensions (Å)", "Longest (Å)", "Charge",
                "Acceptors/Donors", "Folder", "File"]
        self.ligand_table = QTableWidget(0, len(cols))
        self.ligand_table.setHorizontalHeaderLabels(cols)
        self.ligand_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.ligand_table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.ligand_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.ligand_table.setSortingEnabled(True)
        self.ligand_table.setMinimumHeight(170)
        self.ligand_table.horizontalHeader().setStretchLastSection(True)
        self.ligand_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hdr_tips = {3: "Dimensions along the principal axes (largest to smallest), independent of orientation. "
                       "Flexible ligands may extend further.",
                    4: "Largest atom-to-atom distance.",
                    2: "Rotatable bonds.", 5: "Sum of partial charges.",
                    6: "Hydrogen-bond acceptor atoms (OA/NA/SA) / polar hydrogens (HD)."}
        for c, tip in hdr_tips.items():
            self.ligand_table.horizontalHeaderItem(c).setToolTip(tip)
        gl.addWidget(self.ligand_table)
        self.ligand_count = QLabel("No ligands selected.")
        self.ligand_count.setWordWrap(True)
        gl.addWidget(self.ligand_count)
        self.largest_label = QLabel("")
        self.largest_label.setWordWrap(True)
        gl.addWidget(self.largest_label)
        lay.addWidget(g)
        lay.addStretch(1)
        return w

    def _start_dir(self, *candidates):
        for c in candidates:
            if c:
                c = c if os.path.isdir(c) else os.path.dirname(c)
                if c and os.path.isdir(c):
                    return c
        return ""

    # -- receptor
    def _browse_receptor(self):
        p, _ = QFileDialog.getOpenFileName(
            self.tool_window.ui_area, "Select the prepared receptor", self._start_dir(self.receptor_edit.text()),
            "PDBQT files (*.pdbqt);;All files (*)")
        if p:
            self.receptor_edit.setText(p)
            self._receptor_chosen(open_model=True)

    def _check_receptor(self):
        p = self.receptor_edit.text().strip().strip('"')
        if not p:
            self.receptor_status.setText("")
            return False
        msgs = pdbqt_tools.check_receptor_pdbqt(p)
        errs = [m for l, m in msgs if l == "error"]
        warns = [m for l, m in msgs if l == "warning"]
        if errs:
            self.receptor_status.setText("✘ " + " ".join(errs))
            self.receptor_status.setStyleSheet(BAD_STYLE)
            return False
        info = pdbqt_tools.inspect_pdbqt(p)
        text = "✔ %d atoms (%d polar H)." % (info["n_atoms"], info["n_hydrogens"])
        if warns:
            text += " ⚠ " + " ".join(warns)
        self.receptor_status.setText(text)
        self.receptor_status.setStyleSheet(WARN_STYLE if warns else OK_STYLE)
        return True

    def _receptor_chosen(self, open_model=False):
        ok = self._check_receptor()
        path = self.receptor_edit.text().strip().strip('"')
        if ok and open_model and path:
            self._open_receptor(path)

    def _top_models(self):
        return [m for m in self.session.models.list() if len(m.id) == 1]

    def _open_receptor(self, path):
        """Show the receptor in the scene; hide other (unrelated) models so they do not clash."""
        norm = os.path.normcase(os.path.abspath(path))
        if (self._receptor_path == norm and self._receptor_models
                and not any(getattr(m, "deleted", False) for m in self._receptor_models)):
            return
        self._close_receptor_models()
        from .box import DockBoxModel
        keep = [m for m in self._top_models() if isinstance(m, DockBoxModel)]
        hidden = [m for m in self._top_models() if m not in keep and m.display]
        for m in hidden:
            m.display = False
        if hidden:
            self.session.logger.info("Lattice VinaPilotX: hid %d other model(s) (show them again in the Models panel)." % len(hidden))
        before = set(self._top_models())
        try:
            run(self.session, 'open "%s" format pdbqt' % _chimerax_path(path))
        except Exception as err:
            self.session.logger.warning("Lattice VinaPilotX could not open the receptor: %s" % err)
            return
        new = [m for m in self._top_models() if m not in before]
        self._receptor_models = new
        self._receptor_path = norm
        if new:
            try:
                run(self.session, "view #%s" % ",".join(str(m.id[0]) for m in new))
            except Exception:
                pass
            self._preset_box_from_receptor(new)

    def _preset_box_from_receptor(self, models):
        """Pre-fill the Box tab with the receptor's extent (blind-docking starting point); no box is drawn."""
        import numpy as np
        pts = []
        for m in models:
            atoms = getattr(m, "atoms", None)
            if atoms is not None and len(atoms):
                pts.append(np.asarray(atoms.scene_coords, dtype=float))
        if not pts:
            return
        xyz = np.vstack(pts)
        lo, hi = xyz.min(axis=0), xyz.max(axis=0)
        self.box_panel.preset_fields((lo + hi) / 2.0, hi - lo)

    def _close_receptor_models(self):
        for m in self._receptor_models:
            try:
                if not getattr(m, "deleted", False):
                    m.delete()
            except Exception:
                pass
        self._receptor_models = []
        self._receptor_path = None

    # -- ligands
    def _add_ligand_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self.tool_window.ui_area, "Select ligand PDBQT files",
            self._start_dir(self.controller.settings.get("last_ligand_dir"), self.workdir_edit.text()),
            "PDBQT files (*.pdbqt);;All files (*)")
        if files:
            self.controller.settings.set("last_ligand_dir", os.path.dirname(files[0]))
            self._add_ligand_paths(files)

    def _add_ligand_folder(self):
        d = QFileDialog.getExistingDirectory(
            self.tool_window.ui_area, "Select a folder of ligand PDBQT files",
            self._start_dir(self.controller.settings.get("last_ligand_dir"), self.workdir_edit.text()))
        if d:
            self.controller.settings.set("last_ligand_dir", d)
            self._add_ligand_paths([d])

    def _add_ligand_paths(self, paths):
        existing = {r.key for r in self._records.values()}
        QGuiApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            records = pdbqt_tools.scan_ligands(paths, self.recursive_cb.isChecked(), self._used_names,
                                               skip_keys=existing,
                                               progress=lambda i, n: QGuiApplication.processEvents())
        except (OSError, ValueError) as err:
            self.ligand_count.setText("✘ %s" % err)
            return
        finally:
            QGuiApplication.restoreOverrideCursor()
        self.ligand_table.setSortingEnabled(False)
        for rec in records:
            self._add_ligand_row(rec)
        self.ligand_table.setSortingEnabled(True)
        self._update_ligand_summary()

    def _add_ligand_row(self, rec):
        rid = self._next_rid
        self._next_rid += 1
        self._records[rid] = rec
        p = rec.props
        r = self.ligand_table.rowCount()
        self.ligand_table.insertRow(r)
        folder, fname = os.path.split(rec.source)
        fname_label = fname if rec.single else "%s  [#%d]" % (fname, rec.index + 1)
        cells = [
            _NumItem(rec.name), _NumItem(str(p["n_atoms"]), p["n_atoms"]),
            _NumItem("" if p["torsions"] is None else str(p["torsions"]), p["torsions"]),
            _NumItem(_fmt_dims(p["dims"]), p["dims"][0]), _NumItem("%.1f" % p["longest"], p["longest"]),
            _NumItem("%+.2f" % p["charge"], p["charge"]),
            _NumItem("%d / %d" % (p["n_acceptors"], p["n_polar_h"]), p["n_acceptors"]),  # acceptors / donors (HD)
            _NumItem(folder), _NumItem(fname_label),
        ]
        for c, it in enumerate(cells):
            it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
            it.setData(Qt.ItemDataRole.UserRole, rid)
            it.setToolTip(rec.source if c >= 7 else (p["smiles"] if c == 0 and p["smiles"] else ""))
            if rec.problems:
                it.setForeground(QBrush(QColor("#c62828")))
                it.setToolTip("Invalid ligand: " + "; ".join(rec.problems))
            self.ligand_table.setItem(r, c, it)

    def _save_ligand_csv(self):
        if not self._records:
            self.ligand_count.setText("No ligands to save.")
            return
        start = self.workdir_edit.text().strip()
        start = start if os.path.isdir(start) else self._start_dir(self.controller.settings.get("last_ligand_dir"))
        path, _ = QFileDialog.getSaveFileName(self.tool_window.ui_area, "Save ligand table",
                                              os.path.join(start, "ligands.csv"), "CSV (*.csv)")
        if not path:
            return
        ordered = []
        for row in range(self.ligand_table.rowCount()):  # keep the order currently shown
            rec = self._records.get(self.ligand_table.item(row, 0).data(Qt.ItemDataRole.UserRole))
            if rec is not None:
                ordered.append(rec)
        try:
            pdbqt_tools.write_ligand_csv(ordered, path)
        except OSError as err:
            self.ligand_count.setText("✘ Could not save: %s" % err)
            return
        self.ligand_count.setText("Saved %d ligand(s) to %s" % (len(ordered), path))

    def _remove_ligands(self):
        rids = {self.ligand_table.item(i.row(), 0).data(Qt.ItemDataRole.UserRole)
                for i in self.ligand_table.selectedIndexes()}
        rows = sorted({i.row() for i in self.ligand_table.selectedIndexes()}, reverse=True)
        for r in rows:
            self.ligand_table.removeRow(r)
        for rid in rids:
            self._records.pop(rid, None)
        self._update_ligand_summary()

    def _clear_ligands(self):
        self.ligand_table.setRowCount(0)
        self._records.clear()
        self._used_names.clear()
        self._update_ligand_summary()

    def _update_ligand_summary(self):
        recs = list(self._records.values())
        if not recs:
            self.ligand_count.setText("No ligands selected.")
            self.largest_label.setText("")
            return
        bad = [r for r in recs if r.problems]
        files = len({r.source for r in recs})
        text = "%d ligand(s) from %d file(s)." % (len(recs), files)
        if bad:
            text += "  ✘ %d invalid (will be skipped; hover a red row for the reason)." % len(bad)
        self.ligand_count.setText(text)
        good = [r for r in recs if not r.problems]
        largest = max(good, key=lambda r: r.props["longest"], default=None)
        hi = QBrush(QColor(HIGHLIGHT))
        white = QBrush(QColor("#ffffff"))
        none = QBrush()
        font_b = QFont(); font_b.setBold(True)
        for row in range(self.ligand_table.rowCount()):
            rid = self.ligand_table.item(row, 0).data(Qt.ItemDataRole.UserRole)
            is_big = largest is not None and self._records.get(rid) is largest and len(good) > 1
            for c in range(self.ligand_table.columnCount()):
                it = self.ligand_table.item(row, c)
                it.setBackground(hi if is_big else none)
                if is_big:
                    it.setForeground(white)
                elif not self._records[rid].problems:
                    it.setForeground(QBrush())
                it.setFont(font_b if is_big else QFont())
        if largest is not None:
            p = largest.props
            tag = ('<span style="background-color:%s; color:#ffffff;"><b>&nbsp;%s&nbsp;</b></span>')
            self.largest_label.setTextFormat(Qt.TextFormat.RichText)
            self.largest_label.setText(
                "Largest ligand: " + tag % (HIGHLIGHT, "%s &ndash; %s Å" % (largest.name, _fmt_dims(p["dims"]))) +
                "<br>Longest atom-to-atom distance: " + tag % (HIGHLIGHT, "%.1f Å" % p["longest"]) +
                "<br>Suggestion: Ensure the grid box exceeds this in all directions to accommodate full "
                "ligand movement.")
        else:
            self.largest_label.setText("")

    # ---- Run ---------------------------------------------------------
    def _build_run_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox("Vina settings")
        form = QFormLayout(g)
        self.exh_spin = QSpinBox(); self.exh_spin.setRange(1, 1024); self.exh_spin.setValue(DEFAULTS["exhaustiveness"])
        self.exh_spin.setToolTip("Search effort. Default: 8.")
        self.modes_spin = QSpinBox(); self.modes_spin.setRange(1, 100); self.modes_spin.setValue(DEFAULTS["num_modes"])
        self.modes_spin.setToolTip("Max poses per ligand. Default: 9.")
        self.erange_spin = QDoubleSpinBox(); self.erange_spin.setRange(0.1, 100.0)
        self.erange_spin.setValue(DEFAULTS["energy_range"]); self.erange_spin.setSuffix(" kcal/mol")
        self.erange_spin.setToolTip("Energy window for poses. Default: 3 kcal/mol.")
        self.scoring_combo = QComboBox()
        for s in docking.SCORING_FUNCTIONS:
            self.scoring_combo.addItem(s)
        self.scoring_combo.setToolTip("Scoring function. Default: vina.")
        self.cpu_spin = QSpinBox(); self.cpu_spin.setRange(0, 256)
        self.cpu_spin.setSpecialValueText("automatic")
        self.cpu_spin.setToolTip("CPU cores per docking job. Default: automatic (Vina uses all available cores).")
        self.parallel_spin = QSpinBox(); self.parallel_spin.setRange(1, 64); self.parallel_spin.setValue(1)
        self.parallel_spin.setToolTip("Max concurrent docking jobs. Default: 1. Not multi-ligand co-docking.")
        self.seed_spin = QSpinBox(); self.seed_spin.setRange(0, 2147483647)
        self.seed_spin.setSpecialValueText("random")
        self.seed_spin.setToolTip("Default: random. Enter an integer for reproducible runs. "
                                  "Used seeds are saved in logs and tables.")
        form.addRow("Exhaustiveness", self.exh_spin)
        form.addRow("Number of poses", self.modes_spin)
        form.addRow("Energy range", self.erange_spin)
        form.addRow("Scoring function", self.scoring_combo)
        form.addRow("CPU cores per job", self.cpu_spin)
        form.addRow("Parallel jobs", self.parallel_spin)
        form.addRow("Random seed", self.seed_spin)
        self.overwrite_cb = QCheckBox("Re-dock ligands that already have results")
        form.addRow(self.overwrite_cb)
        reset = QPushButton("Reset to defaults")
        reset.clicked.connect(self._reset_defaults)
        form.addRow(reset)
        note = QLabel("Parallel jobs run different docking jobs simultaneously; they do not co-dock ligands.")
        note.setWordWrap(True)
        form.addRow(note)
        lay.addWidget(g)

        self.box_summary = QLabel("")
        self.box_summary.setWordWrap(True)
        lay.addWidget(self.box_summary)

        row = QHBoxLayout()
        self.run_btn = QPushButton("Run docking")
        self.run_btn.setStyleSheet("font-weight: bold; padding: 6px;")
        self.run_btn.clicked.connect(self._run_clicked)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.controller.cancel_docking)
        self.open_dir_btn = QPushButton("Open folder")
        self.open_dir_btn.clicked.connect(self._open_workdir)
        for b in (self.run_btn, self.cancel_btn, self.open_dir_btn):
            row.addWidget(b)
        lay.addLayout(row)

        self.counter_label = QLabel("Ligands completed: 0 / 0")
        f = QFont(); f.setBold(True)
        self.counter_label.setFont(f)
        lay.addWidget(self.counter_label)
        self.slots_box = QVBoxLayout()
        lay.addLayout(self.slots_box)
        self.run_status = QLabel("Ready.")
        self.run_status.setWordWrap(True)
        lay.addWidget(self.run_status)
        self.activity = QPlainTextEdit()
        self.activity.setReadOnly(True)
        self.activity.setMaximumBlockCount(500)
        self.activity.setMinimumHeight(110)
        lay.addWidget(self.activity, 1)
        return w

    def _reset_defaults(self):
        self.exh_spin.setValue(DEFAULTS["exhaustiveness"]); self.modes_spin.setValue(DEFAULTS["num_modes"])
        self.erange_spin.setValue(DEFAULTS["energy_range"]); self.scoring_combo.setCurrentText(DEFAULTS["scoring"])
        self.cpu_spin.setValue(DEFAULTS["cpu"]); self.parallel_spin.setValue(DEFAULTS["parallel_jobs"])
        self.seed_spin.setValue(DEFAULTS["seed"])

    def _params(self):
        return docking.DockingParams(
            exhaustiveness=self.exh_spin.value(), num_modes=self.modes_spin.value(),
            energy_range=self.erange_spin.value(), cpu=self.cpu_spin.value(),
            parallel_jobs=self.parallel_spin.value(), seed=self.seed_spin.value(),
            scoring=self.scoring_combo.currentText())

    def _run_clicked(self):
        box = self.box_panel.current_box()
        if box is None:
            QMessageBox.warning(self.tool_window.ui_area, "Lattice VinaPilotX",
                                "There is no search box yet. Create one on the Box tab first.")
            self.tabs.setCurrentIndex(2)
            return
        center, size = box
        receptor = self.receptor_edit.text().strip().strip('"')
        workdir = self.workdir_edit.text().strip().strip('"')
        records = list(self._records.values())
        ok, result = self.controller.start_docking(
            receptor, records, workdir, center, size, self._params(),
            recursive=self.recursive_cb.isChecked(), overwrite=self.overwrite_cb.isChecked())
        if not ok:
            QMessageBox.warning(self.tool_window.ui_area, "Lattice VinaPilotX – cannot start", str(result))
            return
        self._save_settings()
        self.activity.clear()
        self.tabs.setCurrentIndex(3)

    def _open_workdir(self):
        d = self.workdir_edit.text().strip()
        if d and os.path.isdir(d):
            QDesktopServices.openUrl(QUrl.fromLocalFile(d))

    def _update_box_summary(self):
        box = self.box_panel.current_box()
        if box is None:
            self.box_summary.setText("⚠ No search box yet - create one on the Box tab.")
            self.box_summary.setStyleSheet(WARN_STYLE)
        else:
            c, s = box
            self.box_summary.setText("Search box: center (%.2f, %.2f, %.2f), size %.2f × %.2f × %.2f Å "
                                     "(volume %.0f Å³)." % (*c, *s, docking.box_volume(s)))
            self.box_summary.setStyleSheet("")

    def _tab_changed(self, index):
        if index == 3:
            self._update_box_summary()

    # per-ligand progress bars ("slots")
    def _make_slots(self, n):
        for lab, bar, _ in self._slots:
            lab.setParent(None); bar.setParent(None)
        self._slots = []
        for _ in range(n):
            lab = QLabel("idle")
            bar = QProgressBar(); bar.setRange(0, 100); bar.setValue(0)
            self.slots_box.addWidget(lab); self.slots_box.addWidget(bar)
            self._slots.append([lab, bar, None])

    def _slot_for(self, job, allocate):
        for s in self._slots:
            if s[2] == job.name:
                return s
        if allocate:
            for s in self._slots:
                if s[2] is None:
                    s[2] = job.name
                    return s
        return None

    def _job_progress(self, job):
        if job.status == docking.STATUS_RUNNING:
            s = self._slot_for(job, True)
            if s is not None:
                s[0].setText("%s%s" % (job.name, "" if job.seed is None else "   (seed %s)" % job.seed))
                s[1].setValue(int(job.progress))
        elif job.status in (docking.STATUS_DONE, docking.STATUS_FAILED, docking.STATUS_CANCELLED):
            s = self._slot_for(job, False)
            if s is not None:
                ok = job.status == docking.STATUS_DONE
                s[0].setText("%s - %s" % (job.name, "finished" if ok else job.status))
                s[1].setValue(100 if ok else int(job.progress))
                s[2] = None

    # ---- Results -----------------------------------------------------
    def _build_results_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Ligand", "Best (kcal/mol)", "Poses", "Time (s)", "Seed", "Status"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setMinimumHeight(220)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.cellDoubleClicked.connect(lambda *_: self._open_selected())
        lay.addWidget(self.table, 1)
        self.with_receptor_cb = QCheckBox("Also open the receptor")
        self.with_receptor_cb.setChecked(True)
        lay.addWidget(self.with_receptor_cb)
        row = QHBoxLayout()
        for text, slot in (("Open selected", self._open_selected), ("Open top 5", self._open_top),
                           ("Export CSV…", self._export_csv), ("Load folder…", self._load_folder)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        lay.addLayout(row)
        self.results_note = QLabel("Double-click a row to open its poses.")
        self.results_note.setWordWrap(True)
        lay.addWidget(self.results_note)
        return w

    def _fill_table(self, jobs):
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self._rows = {}
        self.table.setRowCount(len(jobs))
        for r, job in enumerate(jobs):
            items = [QTableWidgetItem(job.name), _NumItem(""), _NumItem("", 0), _NumItem("", 0),
                     _NumItem(""), QTableWidgetItem("")]
            for c, it in enumerate(items):
                it.setFlags(it.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, it)
            self._rows[job.name] = (items, job)
            self._update_row(job)
        self.table.setSortingEnabled(True)

    def _update_row(self, job):
        entry = self._rows.get(job.name)
        if entry is None:
            return
        items, _ = entry
        items[1].setText("" if job.best is None else "%.2f" % job.best); items[1]._value = job.best
        items[2].setText(str(job.n_poses)); items[2]._value = job.n_poses
        items[3].setText("%.1f" % job.seconds if job.seconds else ""); items[3]._value = job.seconds
        items[4].setText("" if job.seed is None else str(job.seed)); items[4]._value = job.seed
        status = job.status
        if job.status == docking.STATUS_RUNNING:
            status = "running %d%%" % job.progress
        elif job.error:
            status = "%s: %s" % (job.status, job.error)
        items[5].setText(status)
        items[5].setToolTip(job.error or "")

    def _selected_jobs(self):
        names = {self.table.item(i.row(), 0).text() for i in self.table.selectedIndexes()}
        return [self._rows[n][1] for n in names if n in self._rows]

    def _open_jobs(self, jobs):
        jobs = [j for j in jobs if os.path.isfile(j.out_path) and j.best is not None]
        if not jobs:
            self.results_note.setText("Nothing to open: select ligands that finished docking.")
            return
        if self.with_receptor_cb.isChecked():
            rec = self._receptor_for_results()
            if rec and not self._receptor_models:
                run(self.session, 'open "%s" format pdbqt' % _chimerax_path(rec))
        for j in sorted(jobs, key=lambda j: j.best):
            run(self.session, 'open "%s" format pdbqt' % _chimerax_path(j.out_path))

    def _receptor_for_results(self):
        c = self.controller
        if c.runner is not None and c.runner.receptor and os.path.isfile(c.runner.receptor):
            return c.runner.receptor
        rdir = os.path.join(self.workdir_edit.text().strip(), "receptor")
        if os.path.isdir(rdir):
            for f in sorted(os.listdir(rdir)):
                if f.lower().endswith(".pdbqt"):
                    return os.path.join(rdir, f)
        return None

    def _open_selected(self):
        self._open_jobs(self._selected_jobs())

    def _open_top(self):
        jobs = [e[1] for e in self._rows.values() if e[1].best is not None]
        jobs.sort(key=lambda j: j.best)
        self._open_jobs(jobs[:5])

    def _export_csv(self):
        if not self._rows:
            return
        path, _ = QFileDialog.getSaveFileName(self.tool_window.ui_area, "Export results",
                                              os.path.join(self.workdir_edit.text() or "", "results_table.csv"),
                                              "CSV (*.csv)")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ligand", "best_affinity_kcal_mol", "n_poses", "runtime_s", "seed", "status", "error"])
            for _, job in sorted(self._rows.values(), key=lambda e: (e[1].best is None, e[1].best or 0)):
                w.writerow([job.name, "" if job.best is None else "%.3f" % job.best, job.n_poses,
                            "%.1f" % job.seconds, "" if job.seed is None else job.seed, job.status, job.error])
        self.results_note.setText("Exported %s" % path)

    def _load_folder(self):
        d = QFileDialog.getExistingDirectory(self.tool_window.ui_area, "Select a Lattice VinaPilotX working directory",
                                             self._start_dir(self.workdir_edit.text()))
        if not d:
            return
        jobs = docking.load_results_dir(d)
        if not jobs:
            self.results_note.setText("No Vina results found in %s/results" % d)
            return
        self.workdir_edit.setText(d)
        self._fill_table(jobs)
        self.results_note.setText("Loaded %d result(s) from %s" % (len(jobs), d))

    # ---------------------------------------------------------- settings
    def _restore_settings(self):
        s = self.controller.settings
        self._restoring = True
        wanted = self.controller.native_pdbqt_wanted
        self.native_cb.setChecked(wanted)
        self._restoring = False
        from . import format_preference
        if wanted and not format_preference.override_active():
            self.controller.set_native_pdbqt(True, quiet=True)  # startup hook did not run or failed
        self.workdir_edit.setText(s.get("last_workdir", "") or "")
        self.receptor_edit.setText(s.get("last_receptor", "") or "")
        if self.receptor_edit.text():
            self._check_receptor()  # checked, but NOT opened until the user chooses it
        p = s.get("params", {}) or {}
        self.exh_spin.setValue(int(p.get("exhaustiveness", DEFAULTS["exhaustiveness"])))
        self.modes_spin.setValue(int(p.get("num_modes", DEFAULTS["num_modes"])))
        self.erange_spin.setValue(float(p.get("energy_range", DEFAULTS["energy_range"])))
        self.cpu_spin.setValue(int(p.get("cpu", 0)))
        self.parallel_spin.setValue(int(p.get("parallel_jobs", 1)))
        self.seed_spin.setValue(int(p.get("seed", 0)))
        i = self.scoring_combo.findText(p.get("scoring", "vina"))
        if i >= 0:
            self.scoring_combo.setCurrentIndex(i)

    def _save_settings(self):
        s = self.controller.settings
        s.data["last_receptor"] = self.receptor_edit.text().strip()
        s.data["last_workdir"] = self.workdir_edit.text().strip()
        from dataclasses import asdict
        s.data["params"] = asdict(self._params())
        s.save()

    # ------------------------------------------------------------ events
    def _on_event(self, event, payload):
        if event == "vina":
            self._refresh_vina_status()
            ok, info = payload
            if not ok and self.controller.vina_exe is None:
                self.vina_status.setText("✘ " + str(info))
                self.vina_status.setStyleSheet(BAD_STYLE)
        elif event == "busy":
            busy = bool(payload)
            if self.controller.downloading:
                self.download_btn.setEnabled(False); self.cancel_dl_btn.setVisible(True)
                self.dl_bar.setVisible(True); self.dl_bar.setRange(0, 0)
            else:
                self.download_btn.setEnabled(True); self.cancel_dl_btn.setVisible(False)
                self.dl_bar.setVisible(False)
            running = self.controller.running or (busy and not self.controller.downloading)
            self.run_btn.setEnabled(not running)
            self.cancel_btn.setEnabled(running)
        elif event == "download_progress":
            done, total = payload
            if total:
                self.dl_bar.setRange(0, 100); self.dl_bar.setValue(int(100 * done / total))
        elif event == "started":
            self._fill_table(payload)
            self._make_slots(min(self.parallel_spin.value(), len(payload)))
            self.counter_label.setText("Ligands completed: 0 / %d" % len(payload))
            self.run_status.setText("Docking %d ligand(s) ..." % len(payload))
        elif event == "job":
            self._update_row(payload)
            self._job_progress(payload)
        elif event == "progress":
            done, total = payload
            self.counter_label.setText("Ligands completed: %d / %d" % (done, total))
        elif event == "log":
            self.activity.appendPlainText(str(payload))
        elif event == "finished":
            self.run_status.setText("Finished. Results are in %s" % payload.workdir)
            self.table.sortItems(1, Qt.SortOrder.AscendingOrder)
            self.tabs.setCurrentIndex(4)
            self.run_btn.setEnabled(True); self.cancel_btn.setEnabled(False)

    # ------------------------------------------------------------ plumbing
    def show_tab(self, index):
        self.tabs.setCurrentIndex(index)

    def delete(self):
        self.controller.remove_listener(self._on_event)
        try:
            self.box_panel.cleanup()
        except RuntimeError:
            pass
        super().delete()

    def take_snapshot(self, session, flags):
        return {}

    @classmethod
    def restore_snapshot(cls, session, data):
        return cls(session, "Lattice VinaPilotX")
