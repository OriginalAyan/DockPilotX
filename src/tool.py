"""The Lattice VinaPilotX tool window: Settings, Prepare, Inputs, Box, Run and Results tabs."""

import csv
import os

from Qt.QtCore import QTimer, Qt, QUrl
from Qt.QtGui import QBrush, QColor, QDesktopServices, QFont, QGuiApplication
from Qt.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QFrame,
    QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QPlainTextEdit,
    QProgressBar, QPushButton, QScrollArea, QSpinBox, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget, QColorDialog,
)
from chimerax.core.tools import ToolInstance
from chimerax.core.commands import run

from . import docking, layout, pdbqt_tools
from .project import Project
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
        self._records = {}            # ligand table: id -> LigandRecord
        self._used_names = set()
        self._next_rid = 0
        self._slots = []              # per-running-ligand progress widgets
        self._receptor_path = None
        self._receptor_models = []
        self._project = None
        self._build_ui()
        self._restore_settings()
        self.controller.add_listener(self._on_event)
        self._refresh_vina_status()
        if self.controller.prep_env is not None and not self.controller.prep_message:
            self.controller.verify_prep()
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
        from .protein_tab import ProteinTab
        from .ligand_tab import LigandTab
        self.protein = ProteinTab(self.session, self.controller, self)
        self.tabs.addTab(_scrolled(self.protein), "2 Protein")
        self.ligands = LigandTab(self.session, self.controller, self)
        self.tabs.addTab(_scrolled(self.ligands), "3 Ligands")
        self.tabs.addTab(_scrolled(self._build_inputs_tab()), "4 Inputs")
        self.box_panel = BoxPanel(self.session)
        self.tabs.addTab(self.box_panel, "5 Box")
        self.tabs.addTab(_scrolled(self._build_run_tab()), "6 Run")
        self.tabs.addTab(_scrolled(self._build_results_tab()), "7 Results")
        self.tabs.currentChanged.connect(self._tab_changed)

    # ---- Settings ----------------------------------------------------
    def _build_settings_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        g = QGroupBox("Setup folder")
        gl = QVBoxLayout(g)
        self.setup_edit = QLineEdit()
        self.setup_edit.setText(self.controller.setup_dir)
        self.setup_edit.setToolTip("Vina, Meeko and PDB2PQR are downloaded and installed inside this folder.")
        self.setup_edit.editingFinished.connect(self._setup_dir_edited)
        row, _ = _row_with_browse(self.setup_edit, "Browse…", self._browse_setup)
        gl.addLayout(row)
        n = QLabel("Everything the tool installs goes into this folder. To uninstall, delete the folder. "
                   "Keep the path short: micromamba fails on very long paths (Windows allows 260 characters). "
                   "A change applies to new installs.")
        n.setWordWrap(True)
        gl.addWidget(n)
        self.setup_warn = QLabel("")
        self.setup_warn.setWordWrap(True)
        self.setup_warn.setStyleSheet(WARN_STYLE)
        gl.addWidget(self.setup_warn)
        lay.addWidget(g)
        self._check_setup_len()

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
        self.workdir_edit.editingFinished.connect(self._workdir_changed)
        row, _ = _row_with_browse(self.workdir_edit, "Browse…", self._browse_workdir)
        gl.addLayout(row)
        n = QLabel("Use a new folder for each docking project. The tool keeps a project file (project.vinapx) and "
                   "folders starting with vinax_ here: do not rename or edit them.")
        n.setWordWrap(True)
        gl.addWidget(n)
        prow = QHBoxLayout()
        self.project_label = QLabel("")
        self.project_label.setWordWrap(True)
        prow.addWidget(self.project_label, 1)
        self.project_open_btn = QPushButton("Open project…")
        self.project_open_btn.setToolTip("Open a project.vinapx file: restores the working directory, receptor, ligands, "
                                         "box and settings.")
        self.project_open_btn.clicked.connect(self._open_project_dialog)
        prow.addWidget(self.project_open_btn)
        gl.addLayout(prow)
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

        g = QGroupBox("Preparation environment (Meeko)")
        gl = QVBoxLayout(g)
        self.env_status = QLabel("")
        self.env_status.setWordWrap(True)
        gl.addWidget(self.env_status)
        self.env_method = QComboBox()
        self.env_method.addItem("Automatic (micromamba, then venv)", "auto")
        self.env_method.addItem("micromamba (conda-forge)", "micromamba")
        self.env_method.addItem("venv + pip (fallback)", "venv")
        self.env_method.setToolTip("Meeko and RDKit are installed in a private environment, not into ChimeraX.")
        gl.addWidget(self.env_method)
        row = QHBoxLayout()
        self.env_setup_btn = QPushButton("Set up Meeko")
        self.env_setup_btn.clicked.connect(lambda: self.controller.setup_prep_env(self.env_method.currentData()))
        self.env_cancel_btn = QPushButton("Cancel")
        self.env_cancel_btn.clicked.connect(self.controller.cancel_prep)
        self.env_cancel_btn.setVisible(False)
        self.env_existing_btn = QPushButton("Use existing Python…")
        self.env_existing_btn.setToolTip("A Python (or environment folder) that already has Meeko and RDKit.")
        self.env_existing_btn.clicked.connect(self._use_existing_python)
        self.env_verify_btn = QPushButton("Verify")
        self.env_verify_btn.clicked.connect(self.controller.verify_prep)
        self.env_remove_btn = QPushButton("Remove")
        self.env_remove_btn.clicked.connect(self._remove_env)
        for b in (self.env_setup_btn, self.env_cancel_btn, self.env_existing_btn, self.env_verify_btn,
                  self.env_remove_btn):
            row.addWidget(b)
        row.addStretch(1)
        gl.addLayout(row)
        self.env_bar = QProgressBar()
        self.env_bar.setRange(0, 0)
        self.env_bar.setVisible(False)
        gl.addWidget(self.env_bar)
        self.env_log = QPlainTextEdit()
        self.env_log.setReadOnly(True)
        self.env_log.setMaximumBlockCount(1000)
        self.env_log.setMinimumHeight(90)
        gl.addWidget(self.env_log)
        lay.addWidget(g)

        g = QGroupBox("Protonation states (PDB2PQR + PROPKA)")
        gl = QVBoxLayout(g)
        self.pq_status = QLabel("")
        self.pq_status.setWordWrap(True)
        gl.addWidget(self.pq_status)
        row = QHBoxLayout()
        self.pq_install_btn = QPushButton("Install PDB2PQR")
        self.pq_install_btn.setToolTip("Installs PDB2PQR (with PROPKA) into the Meeko environment. Set up Meeko first.")
        self.pq_install_btn.clicked.connect(self.controller.install_pdb2pqr)
        row.addWidget(self.pq_install_btn)
        row.addStretch(1)
        gl.addLayout(row)
        lay.addWidget(g)
        self._refresh_env_status()

        from .pose_viewer import DEFAULT_COLORS
        g = QGroupBox("Interaction colours")
        gl = QVBoxLayout(g)
        self._color_btns = {}
        saved = self.controller.settings.get("colors") or {}
        self._colors = {k: saved.get(k, v) for k, v in DEFAULT_COLORS.items()}
        for key, label, tip in (
                ("hbond", "H-bonds", "H-bonds that meet the standard criteria."),
                ("hbond_relaxed", "H-bonds (relaxed criteria)", "H-bonds that only meet the relaxed distance/angle limits."),
                ("clash", "Clashes", "Atom pairs closer than the clash cutoff.")):
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            btn = QPushButton("")
            btn.setFixedWidth(70)
            btn.setToolTip(tip)
            btn.clicked.connect(lambda _=False, k=key: self._pick_color(k))
            self._color_btns[key] = btn
            row.addWidget(btn)
            row.addStretch(1)
            gl.addLayout(row)
        reset = QPushButton("Reset colours")
        reset.clicked.connect(self._reset_colors)
        gl.addWidget(reset)
        self._paint_color_buttons()
        lay.addWidget(g)
        lay.addStretch(1)
        return w

    def _refresh_env_status(self):
        c = self.controller
        if c.prep_ok:
            self.env_status.setText("✔ " + c.prep_message if c.prep_message else "✔ Meeko environment ready.")
            self.env_status.setStyleSheet(OK_STYLE)
        else:
            self.env_status.setText("✘ " + (c.prep_message or "Not set up. Meeko prepares receptors and ligands."))
            self.env_status.setStyleSheet(BAD_STYLE)
        if c.pdb2pqr_ok:
            self.pq_status.setText("✔ PDB2PQR %s ready." % c.prep_info.get("pdb2pqr"))
            self.pq_status.setStyleSheet(OK_STYLE)
        else:
            self.pq_status.setText("✘ Not installed." if c.prep_ok else "✘ Set up Meeko first.")
            self.pq_status.setStyleSheet(BAD_STYLE)

    def _browse_setup(self):
        d = QFileDialog.getExistingDirectory(self.tool_window.ui_area, "Select the setup folder", self.setup_edit.text())
        if d:
            self.setup_edit.setText(d)
            self._setup_dir_edited()

    def _setup_dir_edited(self):
        ok, info = self.controller.set_setup_dir(self.setup_edit.text())
        if ok:
            self.setup_edit.setText(info)
        self._check_setup_len(None if ok else info)

    def _check_setup_len(self, error=None):
        text = error or ""
        if not error and os.name == "nt" and len(self.setup_edit.text()) > 25:
            text = "This path is long. micromamba may fail on Windows; a short folder such as C:\\vpx is safer."
        self.setup_warn.setText(text)
        self.setup_warn.setVisible(bool(text))

    def _prep_busy_changed(self, busy):
        self.env_bar.setVisible(busy)
        self.env_cancel_btn.setVisible(busy)
        for b in (self.env_setup_btn, self.env_existing_btn, self.env_verify_btn, self.env_remove_btn):
            b.setEnabled(not busy)
        self.protein.set_busy(busy)
        self.ligands.set_busy(busy)
        self.pq_install_btn.setEnabled(not busy)

    def _use_existing_python(self):
        p, _ = QFileDialog.getOpenFileName(self.tool_window.ui_area, "Select python (or python.exe)", "",
                                           "All files (*)")
        if p:
            self.controller.use_external_python(p)

    def _remove_env(self):
        if QMessageBox.question(self.tool_window.ui_area, "Remove environment",
                                "Remove the Meeko environment that Lattice VinaPilotX created?") == QMessageBox.StandardButton.Yes:
            self.controller.remove_prep_env()

    def _paint_color_buttons(self):
        for k, btn in self._color_btns.items():
            btn.setStyleSheet("background-color: %s; border: 1px solid #666;" % self._colors[k])

    def _apply_colors(self):
        self.controller.settings.set("colors", dict(self._colors))
        self._paint_color_buttons()
        if hasattr(self, "results"):
            self.results.viewer.set_colors(self._colors)

    def _pick_color(self, key):
        c = QColorDialog.getColor(QColor(self._colors[key]), self.tool_window.ui_area, "Choose colour")
        if c.isValid():
            self._colors[key] = c.name()
            self._apply_colors()

    def _reset_colors(self):
        from .pose_viewer import DEFAULT_COLORS
        self._colors = dict(DEFAULT_COLORS)
        self._apply_colors()

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
            self._workdir_changed()

    # ---- project file ------------------------------------------------
    def get_project(self):
        """The Project of the working directory (None if no folder is set)."""
        wd = self.workdir_edit.text().strip().strip('"')
        if not wd:
            return None
        wd = os.path.abspath(wd)
        if self._project is None or self._project.workdir != wd:
            self._project = Project(wd)
            self._project.load()
        return self._project

    def _workdir_changed(self):
        proj = self.get_project()
        if proj is None:
            self.project_label.setText("")
            return
        self.protein.load_project()
        if proj.loaded:
            fresh = not self._records and not self.receptor_edit.text().strip()
            self.project_label.setText("✔ Project file found (%d receptor revision(s))." % len(proj.revisions))
            self.project_label.setStyleSheet(OK_STYLE)
            if fresh:
                self.restore_project(proj)
        else:
            self.project_label.setText("New project folder. A project.vinapx file will be saved here as you work.")
            self.project_label.setStyleSheet("")

    def save_project(self, last_run=None):
        """Write the inputs, box and parameters into project.vinapx (never raises)."""
        try:
            proj = self.get_project()
            if proj is None:
                return False
            from dataclasses import asdict
            box = self.box_panel.current_box()
            proj.set_inputs(ligand_paths=[r.path for r in self._records.values()], box=box,
                            params=asdict(self._params()), last_run=last_run)
            return proj.save(self.controller_version())
        except Exception:
            return False

    @staticmethod
    def controller_version():
        return docking.TOOL_VERSION

    def restore_project(self, proj):
        """Bring the tool back to the state saved in *proj*."""
        self.workdir_edit.setText(proj.workdir)
        self._project = proj
        self.protein.load_project()
        pdbqt = proj.pdbqt()
        if pdbqt:
            self.receptor_edit.setText(pdbqt)
            self._check_receptor()
        paths = [p for p in proj.ligand_paths() if os.path.exists(p)]
        if paths:
            self._add_ligand_paths(paths)
        box = proj.data.get("box")
        if box:
            try:
                self.box_panel.set_box_from_cli(tuple(box["center"]), tuple(box["size"]), self.box_panel.spacing.value())
            except Exception:
                pass
        self._apply_params(proj.data.get("params") or {})
        self.project_label.setText("✔ Restored the project in %s (%d receptor revision(s), %d ligand(s))." % (
            proj.workdir, len(proj.revisions), len(paths)))
        self.project_label.setStyleSheet(OK_STYLE)
        self.session.logger.info("Lattice VinaPilotX: restored the project in %s" % proj.workdir)

    def open_project(self, path):
        """Open a project.vinapx file (or a folder that contains one)."""
        proj = Project.open(path)
        self.restore_project(proj)
        return proj

    def _open_project_dialog(self):
        p, _ = QFileDialog.getOpenFileName(self.tool_window.ui_area, "Open a project", self._start_dir(self.workdir_edit.text()),
                                           "Lattice VinaPilotX project (*.vinapx);;All files (*)")
        if p:
            try:
                self.open_project(p)
            except Exception as err:
                self.project_label.setText("✘ %s" % err)
                self.project_label.setStyleSheet(BAD_STYLE)

    def drop_prepared_receptor(self, reason):
        """Clear the receptor field if it holds a receptor prepared by this tool, so a new protein is never mixed
        up with an older one."""
        cur = self.receptor_edit.text().strip().strip('"')
        proj = self.get_project()
        if not cur or proj is None:
            return False
        inside = os.path.normcase(os.path.abspath(cur)).startswith(
            os.path.normcase(os.path.abspath(layout.path(proj.workdir, "prepared"))))
        if not inside:
            return False
        self.receptor_edit.setText("")
        self.receptor_status.setText("")
        self._close_receptor_models()
        self.session.logger.info("Lattice VinaPilotX: receptor cleared on the Inputs tab (%s)." % reason)
        return True

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

    def _schedule_project_save(self):
        """Save the project a moment after the last change (several changes in a row give one write)."""
        if getattr(self, "_save_pending", False):
            return
        self._save_pending = True

        def go():
            self._save_pending = False
            self.save_project()
        QTimer.singleShot(1500, go)

    def _receptor_chosen(self, open_model=False):
        self._schedule_project_save()
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
        self._schedule_project_save()
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
            self.tabs.setCurrentIndex(4)
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
        self.save_project(last_run={"workdir": workdir, "receptor": receptor, "ligands": len(records)})
        self.activity.clear()
        self.tabs.setCurrentIndex(5)

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
        if index == 5:
            self._update_box_summary()
        elif index == 1:
            self.protein.refresh_models()

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
        from .results_tab import ResultsTab
        self.results = ResultsTab(self.session, self.controller, self)
        return self.results

    # PoseViewer host interface: the receptor is opened once and reused
    @property
    def receptor_path_hint(self):
        return self._receptor_path

    def get_receptor_models(self):
        return list(self._receptor_models)

    def set_receptor_models(self, models, path):
        self._receptor_models = list(models)
        self._receptor_path = os.path.normcase(os.path.abspath(path))

    def close_receptor(self):
        self._close_receptor_models()

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
        self._apply_params(s.get("params", {}) or {})

    def _apply_params(self, p):
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
            runner = self.controller.runner
            self.results.set_jobs(payload, runner.workdir, runner.receptor)
            self._make_slots(min(self.parallel_spin.value(), len(payload)))
            self.counter_label.setText("Ligands completed: 0 / %d" % len(payload))
            self.run_status.setText("Docking %d ligand(s) ..." % len(payload))
        elif event == "job":
            self.results.update_job(payload)
            self._job_progress(payload)
        elif event == "progress":
            done, total = payload
            self.counter_label.setText("Ligands completed: %d / %d" % (done, total))
        elif event == "log":
            self.activity.appendPlainText(str(payload))
        elif event == "setup_log":
            self.env_log.appendPlainText(str(payload))
        elif event == "protein_log":
            self.protein.log_line(payload)
        elif event == "ligand_log":
            self.ligands.log_line(payload)
        elif event == "prep_busy":
            self._prep_busy_changed(bool(payload))
        elif event == "prep_env":
            self.protein.refresh_env()
            self.ligands.refresh_env()
            self._refresh_env_status()
        elif event == "prep_progress":
            self.ligands.progress(payload)
        elif event == "prep_receptor_done":
            self.protein.receptor_done(payload)
        elif event == "pdb2pqr_done":
            self.protein.pdb2pqr_done(payload)
        elif event == "prep_ligands_done":
            self.ligands.ligands_done(payload)
        elif event == "finished":
            self.run_status.setText("Finished. Results are in %s" % payload.workdir)
            self.results.sort_by_best()
            self.tabs.setCurrentIndex(6)
            self.run_btn.setEnabled(True); self.cancel_btn.setEnabled(False)

    # ------------------------------------------------------------ plumbing
    def show_tab(self, index):
        self.tabs.setCurrentIndex(index)

    def delete(self):
        self.save_project()
        self.controller.remove_listener(self._on_event)
        try:
            self.results.close()
        except Exception:
            pass
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
