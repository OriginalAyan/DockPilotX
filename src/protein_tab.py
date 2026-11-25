"""Protein tab: clean a PDB/mmCIF receptor with Dock Prep, set hydrogens, minimize, and convert it with Meeko.

Every step writes a checkpoint file into the receptor's *revision* folder and is recorded in the project file,
so a project can be closed and continued later, and a new start never mixes with an older protein.
"""

import os
import shutil

from Qt.QtCore import QTimer
from Qt.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
    QProgressBar, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)

import re

from . import layout, prep_undo
from .pdbqt_tools import safe_name
from .prep_widgets import BAD_STYLE, OK_STYLE, WARN_STYLE, EnvBanner, PkaTable
from .receptor_clean import SIDE_CHAIN_CHOICES, STANDARDIZE_NOTE

MIN_TEXT = "3. Minimize (optional)"
POLAR_MODES = [
    ("Keep as in the file (add none)", "keep"),
    ("ChimeraX (residue names, else H-bond method)", "chimerax"),
    ("PDB2PQR + PROPKA at pH", "pdb2pqr"),
]


def stem(name):
    """Model name without a file extension ('cleaned.pdb' -> 'cleaned')."""
    return re.sub(r"\.(pdb|ent|cif|mmcif|pdbqt|mol2|sdf)$", "", str(name or ""), flags=re.I)


def _small_note(text):
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setStyleSheet("color: gray;")
    return lbl


class ProteinTab(QWidget):
    def __init__(self, session, controller, host):
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host          # the tool: workdir_edit, receptor_edit, _receptor_chosen, show_tab, get_project
        self._model_ids = []
        self._prep_model = None       # the cleaned copy in the scene
        self._prep_name = None
        self._rev = None              # active revision (dict of the project file)
        self._orig_model = None       # open model the copy was made from (hidden while the copy is shown)
        self._min_watch = None
        self._h_pre = None            # checkpoint saved before the running hydrogen step
        self._h_ends = []
        self._h_in_pdb = None
        self._swap_state = None
        self._convert_rev = None
        self._handlers = []
        self._refresh_pending = False
        self._build_ui()
        self.refresh_env()
        self._hook_models()
        self.refresh_models()

    # keep the 'Use open model' list current while the window stays open
    def _hook_models(self):
        try:
            from chimerax.core.models import ADD_MODELS, REMOVE_MODELS
            trig = self.session.triggers
            for name in (ADD_MODELS, REMOVE_MODELS):
                self._handlers.append(trig.add_handler(name, self._models_changed))
        except Exception:
            self._handlers = []

    def _models_changed(self, *_):
        if self._refresh_pending:
            return
        self._refresh_pending = True
        QTimer.singleShot(0, self._deferred_refresh)

    def _deferred_refresh(self):
        self._refresh_pending = False
        try:
            self.refresh_models()
        except RuntimeError:      # the window is already gone
            pass

    def cleanup(self):
        for h in self._handlers:
            try:
                self.session.triggers.remove_handler(h)
            except Exception:
                pass
        self._handlers = []

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        lay = QVBoxLayout(self)
        self.banner = EnvBanner(self.controller, self.host)
        lay.addWidget(self.banner)

        g = QGroupBox("Saved receptor preparations")
        gl = QVBoxLayout(g)
        self.rev_combo = QComboBox()
        self.rev_combo.setToolTip("Each start from a protein is saved as a revision in the working directory.")
        gl.addWidget(self.rev_combo)
        row = QHBoxLayout()
        self.rev_continue_btn = QPushButton("Continue from latest step")
        self.rev_continue_btn.setToolTip("Opens the last saved structure of this revision so you can carry on.")
        self.rev_continue_btn.clicked.connect(self._continue_revision)
        row.addWidget(self.rev_continue_btn)
        row.addStretch(1)
        gl.addLayout(row)
        self.rev_label = QLabel("")
        self.rev_label.setWordWrap(True)
        gl.addWidget(self.rev_label)
        lay.addWidget(g)

        g = QGroupBox("Receptor (PDB / mmCIF)")
        gl = QVBoxLayout(g)
        row = QHBoxLayout()
        self.rec_edit = QLineEdit()
        self.rec_edit.setPlaceholderText("Protein file, e.g. 1abc.pdb or 1abc.cif")
        row.addWidget(self.rec_edit, 1)
        b = QPushButton("Browse…")
        b.clicked.connect(self._browse_receptor)
        row.addWidget(b)
        gl.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel("Use open model:"))
        self.model_combo = QComboBox()
        self.model_combo.setToolTip("Prepare a copy of a model that is already open. The model itself is not changed "
                                    "(it is hidden while the copy is shown).")
        row.addWidget(self.model_combo, 1)
        self.model_refresh_btn = QPushButton("⟳")
        self.model_refresh_btn.setFixedWidth(32)
        self.model_refresh_btn.setToolTip("Refresh the list of open models.")
        self.model_refresh_btn.clicked.connect(lambda: self.refresh_models())
        row.addWidget(self.model_refresh_btn)
        gl.addLayout(row)
        self.chains_edit = QLineEdit()
        self.chains_edit.setPlaceholderText("Keep only chains, e.g. A,B (blank = all)")
        gl.addWidget(self.chains_edit)
        lay.addWidget(g)

        # ---- 1 Clean (Dock Prep wording)
        g = QGroupBox("1. Clean (Dock Prep)")
        gl = QVBoxLayout(g)
        self.water_cb = QCheckBox("Delete solvent")
        self.ions_cb = QCheckBox("Delete non-complexed ions")
        self.alt_cb = QCheckBox("Delete non-current alternate locations")
        self.bb_cb = QCheckBox("Delete residues with incomplete backbones")
        self.lig_cb = QCheckBox("Delete other ligands / cofactors")
        self.std_cb = QCheckBox("\"Standardize\" certain residue types:")
        tips = {
            self.water_cb: "Dock Prep option.",
            self.ions_cb: "Dock Prep option. Ions not bound to the ligands are removed.",
            self.alt_cb: "Dock Prep option. Keeps the alternate location in use (highest occupancy).",
            self.bb_cb: "Dock Prep option (off by default).",
            self.lig_cb: "Not part of Dock Prep. Removes ligands, cofactors and other non-polymer groups.",
            self.std_cb: "Dock Prep option.",
        }
        for cb in (self.water_cb, self.ions_cb, self.alt_cb, self.bb_cb, self.lig_cb, self.std_cb):
            cb.setChecked(cb is not self.bb_cb)
            cb.setToolTip(tips[cb])
            gl.addWidget(cb)
            if cb is self.std_cb:
                note = QLabel("\n".join("   • " + t for t in STANDARDIZE_NOTE))
                note.setStyleSheet("color: gray;")
                gl.addWidget(note)
        row = QHBoxLayout()
        self.sc_cb = QCheckBox("Incomplete side chains:")
        self.sc_cb.setChecked(True)
        self.sc_cb.setToolTip("Dock Prep option. Unchecked: incomplete side chains are left as they are.")
        row.addWidget(self.sc_cb)
        self.sc_combo = QComboBox()
        for label, value in SIDE_CHAIN_CHOICES:
            self.sc_combo.addItem(label, value)
        self.sc_cb.toggled.connect(self.sc_combo.setEnabled)
        row.addWidget(self.sc_combo, 1)
        gl.addLayout(row)
        self.rec_btn = QPushButton("1. Clean receptor")
        self.rec_btn.setToolTip("Runs Dock Prep on a copy and keeps it in the scene as <name>_prep. Adds no hydrogens. "
                                "Starts a new saved revision.")
        self.rec_btn.clicked.connect(self._clean_receptor)
        gl.addWidget(self.rec_btn)
        lay.addWidget(g)

        # ---- 2 Hydrogens
        g = QGroupBox("2. Hydrogens")
        gl = QVBoxLayout(g)
        gl.addWidget(QLabel("Polar hydrogens (protonation states)"))
        self.polar_combo = QComboBox()
        for label, value in POLAR_MODES:
            self.polar_combo.addItem(label, value)
        self.polar_combo.setToolTip(
            "Keep: add no polar hydrogens; the ones in the file stay as they are.\n"
            "ChimeraX: residue names are used (HID/HIE/HIP, ASH, GLH, LYN, CYM); plain HIS is decided by H-bonds;\n"
            "   other residues get the usual pH ~7 state.\n"
            "PDB2PQR: PROPKA predicts pKa values and the states follow the pH you enter. Existing hydrogens are replaced.\n"
            "   Chain ends follow PROPKA when they are real ends and are neutral when the chain was cut.")
        row = QHBoxLayout()
        row.addWidget(self.polar_combo, 1)
        self.ph_spin = QDoubleSpinBox()
        self.ph_spin.setRange(0.0, 14.0)
        self.ph_spin.setDecimals(2)
        self.ph_spin.setSingleStep(0.1)
        self.ph_spin.setValue(7.0)
        self.ph_spin.setPrefix("pH ")
        self.ph_spin.setToolTip("pH used to set the protonation states (PDB2PQR only). Default 7.0.")
        row.addWidget(self.ph_spin)
        gl.addLayout(row)
        self.polar_combo.currentIndexChanged.connect(self._polar_changed)
        self.nonpolar_cb = QCheckBox("Add non-polar hydrogens (on carbon)")
        self.nonpolar_cb.setToolTip("Separate from the polar ones. Docking does not need them; minimization does.\n"
                                    "Only hydrogens on carbon are kept, so polar hydrogens you removed stay removed.")
        gl.addWidget(self.nonpolar_cb)
        row = QHBoxLayout()
        self.h_btn = QPushButton("2. Apply hydrogens")
        self.h_btn.setToolTip("Applies the choices above to the cleaned copy (it is cleaned first if needed).")
        self.h_btn.clicked.connect(self._apply_hydrogens)
        row.addWidget(self.h_btn)
        self.h_undo_btn = QPushButton("Undo")
        self.h_undo_btn.setToolTip("Go back to the structure before the last hydrogen step. (Also works with ChimeraX's undo.)")
        self.h_undo_btn.clicked.connect(self._undo_hydrogens)
        self.h_undo_btn.setEnabled(False)
        row.addWidget(self.h_undo_btn)
        self.h_redo_btn = QPushButton("Redo")
        self.h_redo_btn.clicked.connect(self._redo_hydrogens)
        self.h_redo_btn.setEnabled(False)
        row.addWidget(self.h_redo_btn)
        gl.addLayout(row)
        self.h_result = QLabel("")
        self.h_result.setWordWrap(True)
        gl.addWidget(self.h_result)
        self.pka_table = PkaTable()
        gl.addWidget(self.pka_table)
        lay.addWidget(g)
        self._polar_changed()

        # ---- 3 Minimize, 4 Convert
        g = QGroupBox("3 / 4. Minimize and convert")
        gl = QVBoxLayout(g)
        row = QHBoxLayout()
        self.min_btn = QPushButton(MIN_TEXT)
        self.min_btn.setToolTip("ChimeraX energy minimization of the cleaned copy. Moves all atoms. "
                                "Charges are assigned first (addcharge). ChimeraX reports energies every 100 steps.")
        self.min_btn.clicked.connect(self._minimize)
        row.addWidget(self.min_btn)
        self.min_steps = QSpinBox()
        self.min_steps.setRange(1, 100000)
        self.min_steps.setValue(100)
        self.min_steps.setSuffix(" steps")
        row.addWidget(self.min_steps)
        self.min_conv = QCheckBox("until convergence")
        self.min_conv.setToolTip("Not recommended: can take very long for a whole protein.")
        row.addWidget(self.min_conv)
        gl.addLayout(row)
        gl.addWidget(_small_note("Minimization needs the non-polar hydrogens too (step 2). They are not added automatically."))
        self.bad_cb = QCheckBox("Skip residues Meeko cannot handle")
        self.bad_cb.setToolTip("Residues without a matching template are dropped, with a warning, instead of stopping.")
        self.bad_cb.setChecked(True)
        gl.addWidget(self.bad_cb)
        self.conv_btn = QPushButton("4. Convert with Meeko")
        self.conv_btn.setToolTip("Saves the current <name>_prep model and converts it to PDBQT. "
                                 "If step 1 was not run, it is run first.")
        self.conv_btn.clicked.connect(self._convert_receptor)
        gl.addWidget(self.conv_btn)
        self.approve_btn = QPushButton("Approve for docking")
        self.approve_btn.setToolTip("Only approved receptors are listed on the Inputs tab. Inspect the Meeko result first.")
        self.approve_btn.setEnabled(False)
        self.approve_btn.clicked.connect(self._toggle_approve)
        gl.addWidget(self.approve_btn)
        self.rec_result = QLabel("")
        self.rec_result.setWordWrap(True)
        gl.addWidget(self.rec_result)
        lay.addWidget(g)

        row = QHBoxLayout()
        self.bar = QProgressBar()
        self.bar.setVisible(False)
        row.addWidget(self.bar, 1)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self.controller.cancel_prep)
        row.addWidget(self.cancel_btn)
        lay.addLayout(row)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(2000)
        self.log.setMinimumHeight(110)
        lay.addWidget(self.log)
        lay.addStretch(1)

    # --------------------------------------------------------------- state
    def refresh_env(self):
        self.banner.refresh()
        self.set_busy(self.controller.prep_busy)

    def set_busy(self, busy):
        ok = self.controller.prep_ok and not busy
        self.conv_btn.setEnabled(ok)
        self.rec_btn.setEnabled(not busy)
        self.h_btn.setEnabled(not busy)
        self.min_btn.setEnabled(True)
        self.cancel_btn.setVisible(busy)
        self.bar.setVisible(busy)
        if busy:
            self.bar.setRange(0, 0)

    def log_line(self, text):
        self.log.appendPlainText(str(text))

    def refresh_models(self):
        """List the open atomic models for 'Use open model' and keep the current choice if it is still open."""
        keep = None
        i = self.model_combo.currentIndex()
        if 0 < i < len(self._model_ids):
            keep = self._model_ids[i]
        self.model_combo.blockSignals(True)
        self.model_combo.clear()
        self.model_combo.addItem("(none: use the file above)")
        self._model_ids = [None]
        try:
            from chimerax.atomic import AtomicStructure
            for m in self.session.models.list(type=AtomicStructure):
                if getattr(m, "deleted", False):
                    continue
                self.model_combo.addItem("#%s  %s" % (".".join(str(i) for i in m.id), m.name))
                self._model_ids.append(m)
        except Exception:
            pass
        if keep is not None and keep in self._model_ids:
            self.model_combo.setCurrentIndex(self._model_ids.index(keep))
        self.model_combo.blockSignals(False)

    def _polar_changed(self, *_):
        self.ph_spin.setEnabled(self.polar_combo.currentData() == "pdb2pqr")

    # ------------------------------------------------------------ project / revisions
    def _project(self):
        return self.host.get_project()

    def load_project(self):
        """Fill the saved-preparations list from the project file of the working directory."""
        proj = self._project()
        self.rev_combo.clear()
        revs = proj.revisions if proj else []
        for r in revs:
            self.rev_combo.addItem(proj.summary(r), r["id"])
        if proj and proj.data["receptor"]["active"]:
            i = self.rev_combo.findData(proj.data["receptor"]["active"])
            if i >= 0:
                self.rev_combo.setCurrentIndex(i)
        self.rev_label.setText("" if revs else "None yet. Step 1 starts the first one.")
        self.rev_continue_btn.setEnabled(bool(revs))
        self._update_approve()

    def _selected_revision(self):
        proj = self._project()
        return (proj, proj.revision(self.rev_combo.currentData())) if proj and self.rev_combo.currentData() else (proj, None)

    def _continue_revision(self):
        from .receptor_clean import open_prep_model
        proj, rev = self._selected_revision()
        if rev is None:
            return
        stage, path = proj.latest_stage(rev)
        if path is None:
            self._fail("This revision has no saved structure. Run step 1 again.")
            return
        try:
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            self._prep_model = open_prep_model(self.session, path, rev["name"])
        except Exception as err:
            self._fail("Could not open %s: %s" % (path, err))
            return
        self._rev, self._prep_name = rev, rev["name"]
        self._update_approve()
        proj.data["receptor"]["active"] = rev["id"]
        proj.save()
        self.pka_table.setVisible(False)
        self.rec_result.setText("✔ Continuing %s from its last saved step (%s). Run the next step." % (rev["id"], stage))
        self.rec_result.setStyleSheet(OK_STYLE)
        self.log_line("Opened %s (%s step) as %s_prep." % (os.path.basename(path), stage, rev["name"]))

    def _update_approve(self):
        """Approve button: enabled when the active revision has a converted PDBQT."""
        proj, rev = self._project(), self._rev
        ok = bool(proj and rev and proj.pdbqt(rev))
        self.approve_btn.setEnabled(ok)
        self.approve_btn.setText("Remove approval" if ok and rev.get("approved") else "Approve for docking")

    def _toggle_approve(self):
        proj, rev = self._project(), self._rev
        if not (proj and rev and proj.pdbqt(rev)):
            return
        now = proj.approve(rev, not rev.get("approved"))
        proj.save()
        self._update_approve()
        self.host.refresh_prepared_list()
        if now:
            self.rec_result.setText("✔ %s approved: it is now in the receptor list on the Inputs tab." % os.path.basename(proj.pdbqt(rev)))
            self.rec_result.setStyleSheet(OK_STYLE)
        else:
            self.rec_result.setText("Approval removed. It is no longer in the receptor list on the Inputs tab.")
            self.rec_result.setStyleSheet("")

    def _checkpoint(self, stage, **info):
        """Save the current copy as the checkpoint of *stage* and record it in the project file."""
        from .receptor_clean import save_model_pdb
        proj = self._project()
        if proj is None or self._rev is None or not self._prep_alive():
            return None
        path = proj.stage_file(self._rev, stage)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        save_model_pdb(self.session, self._prep_model, path)
        proj.set_stage(self._rev, stage, path, **info)
        proj.save()
        self.load_project()
        return path

    # ------------------------------------------------------------ helpers
    def _workdir(self):
        w = self.host.workdir_edit.text().strip().strip('"')
        if not w:
            self.rec_result.setText("Choose a working directory on the Settings tab first.")
            self.rec_result.setStyleSheet(BAD_STYLE)
            return None
        return os.path.abspath(w)

    def _browse_receptor(self):
        p, _ = QFileDialog.getOpenFileName(self, "Select a protein structure", self.rec_edit.text(),
                                           "Structures (*.pdb *.ent *.cif *.mmcif);;All files (*)")
        if p:
            self.rec_edit.setText(p)
            self.model_combo.setCurrentIndex(0)

    def _source(self, wd):
        """(source file, name, open model or None) from the file field or the selected open model."""
        idx = self.model_combo.currentIndex()
        if 0 < idx < len(self._model_ids):
            from chimerax.core.commands import run
            m = self._model_ids[idx]
            name = safe_name(stem(m.name)) or "receptor"
            src_dir = os.path.join(layout.make(wd, "prepared"), "sources")
            os.makedirs(src_dir, exist_ok=True)
            source = os.path.join(src_dir, name + "_source.pdb")
            spec = "#" + ".".join(str(i) for i in m.id)
            run(self.session, 'save "%s" %s relModel %s' % (source.replace("\\", "/"), spec, spec))
            return source, name, m
        source = self.rec_edit.text().strip().strip('"')
        if not source or not os.path.isfile(source):
            return None
        return source, safe_name(os.path.splitext(os.path.basename(source))[0]) or "receptor", None

    def _fail(self, text, label=None):
        (label or self.rec_result).setText(text)
        (label or self.rec_result).setStyleSheet(BAD_STYLE)

    def _prep_alive(self):
        m = self._prep_model
        return m is not None and not getattr(m, "deleted", False)

    # ------------------------------------------------------- 1 clean
    def _clean_receptor(self):
        """Step 1: Dock Prep on a copy; the copy stays in the scene and starts a new revision."""
        from .receptor_clean import clean_receptor
        wd = self._workdir()
        if wd is None:
            return False
        try:
            src = self._source(wd)
            if src is None:
                self._fail("Choose a protein file or an open model.")
                return False
            source, name, orig = src
            proj = self._project()
            options = {"delete_solvent": self.water_cb.isChecked(), "delete_ions": self.ions_cb.isChecked(),
                       "delete_altlocs": self.alt_cb.isChecked(), "delete_incomplete_backbone": self.bb_cb.isChecked(),
                       "delete_other_ligands": self.lig_cb.isChecked(), "standardize": self.std_cb.isChecked(),
                       "side_chains": self.sc_combo.currentData() if self.sc_cb.isChecked() else "",
                       "chains": self.chains_edit.text().strip()}
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            model, summary = clean_receptor(
                self.session, source, name, chains=self.chains_edit.text(),
                del_water=options["delete_solvent"], del_ions=options["delete_ions"],
                del_ligands=options["delete_other_ligands"], del_altlocs=options["delete_altlocs"],
                del_incomplete=options["delete_incomplete_backbone"], standardize=options["standardize"],
                side_chains=options["side_chains"], log=self.log_line)
        except Exception as err:
            self._fail("Cleaning failed: %s" % err)
            return False
        self._prep_model, self._prep_name = model, name
        self._rev = proj.new_revision(name, proj.rel(source) if orig is None else "open model " + str(orig.name), options)
        self._checkpoint("clean")
        if orig is not None and not getattr(orig, "deleted", False):
            try:
                orig.display = False
                self.log_line("Hid the original model #%s (%s); the cleaned copy is shown instead."
                              % (".".join(str(i) for i in orig.id), orig.name))
            except Exception:
                pass
        try:
            model.display = True
        except Exception:
            pass
        self._orig_model = orig
        self.pka_table.setVisible(False)
        self._set_undo_state(None)
        self.log_line(summary)
        self.rec_result.setText("✔ %s Saved as revision %s. Step 2 sets hydrogens (optional)." % (summary, self._rev["id"]))
        self.rec_result.setStyleSheet(OK_STYLE)
        return True

    # ------------------------------------------------------- 2 hydrogens
    def _apply_hydrogens(self):
        from .receptor_clean import save_model_pdb
        mode = self.polar_combo.currentData()
        if mode == "pdb2pqr" and not self.controller.pdb2pqr_ok:
            self._fail("PDB2PQR is not installed. Install it on the Settings tab (below Meeko).", self.h_result)
            return
        if not self._prep_alive() and not self._clean_receptor():
            return
        wd = self._workdir()
        if wd is None:
            return
        self.h_result.setText("")
        # checkpoint of the structure before this step, for undo
        proj = self._project()
        self._h_pre = os.path.join(proj.rev_dir(self._rev), "_before_hydrogens.pdb")
        try:
            save_model_pdb(self.session, self._prep_model, self._h_pre)
        except Exception as err:
            self._fail("Saving failed: %s" % err, self.h_result)
            return
        if mode == "pdb2pqr":
            self._start_pdb2pqr()
            return
        try:
            notes = []
            if mode == "chimerax":
                from .hydrogens import add_hydrogens
                d = add_hydrogens(self.session, self._prep_model, "polar", self.log_line)
                notes.append("polar hydrogens: %d kept (ChimeraX placed %d in total, %d on carbon removed)"
                             % (d["kept"], d["placed"], d["removed"]))
            self._finish_hydrogens(notes, mode)
        except Exception as err:
            self._fail("Adding hydrogens failed: %s" % err, self.h_result)

    def _start_pdb2pqr(self):
        from .receptor_clean import save_model_pdb
        from .hydrogens import chain_ends
        proj = self._project()
        out_dir = os.path.join(proj.rev_dir(self._rev), "pdb2pqr")
        os.makedirs(out_dir, exist_ok=True)
        in_pdb = os.path.join(out_dir, self._prep_name + "_before_pdb2pqr.pdb")
        try:
            save_model_pdb(self.session, self._prep_model, in_pdb)
            self._h_ends = chain_ends(self._prep_model)       # real or cut ends, from the structure with SEQRES
        except Exception as err:
            self._fail("Saving failed: %s" % err, self.h_result)
            return
        n_art = sum(1 for e in self._h_ends if not e["real"])
        self.log_line("Chain ends: %d real, %d cut (neutral)." % (len(self._h_ends) - n_art, n_art))
        self._h_in_pdb = in_pdb
        self.h_result.setText("Running PDB2PQR + PROPKA ...")
        self.h_result.setStyleSheet("")
        self.controller.run_pdb2pqr_async(in_pdb, out_dir, self._prep_name, self.ph_spin.value())

    def _finish_hydrogens(self, notes, mode):
        """Non-polar step (if wanted), checkpoint, undo registration and the summary line."""
        from .hydrogens import add_hydrogens, count_hydrogens
        if self.nonpolar_cb.isChecked():
            d = add_hydrogens(self.session, self._prep_model, "nonpolar", self.log_line)
            notes.append("non-polar hydrogens: %d kept (ChimeraX placed %d in total, %d on N/O/S removed)"
                         % (d["kept"], d["placed"], d["removed"]))
        pol, nonpol = count_hydrogens(self._prep_model)
        notes.append("now %d polar and %d non-polar hydrogens" % (pol, nonpol))
        after = self._checkpoint("hydrogens", polar_mode=mode, ph=self.ph_spin.value() if mode == "pdb2pqr" else None,
                                 nonpolar=self.nonpolar_cb.isChecked())
        if self._h_pre and after:
            act = prep_undo.ModelSwapUndo("receptor hydrogens", self._h_pre, after, self._swap_model)
            prep_undo.register(self.session, act)
            self._set_undo_state("after")
        text = "✔ " + "; ".join(notes) + "."
        self.log_line(text)
        self.h_result.setText(text)
        self.h_result.setStyleSheet(OK_STYLE)

    def pdb2pqr_done(self, payload):
        from . import pdb2pqr_tools as pq
        from .hydrogens import fix_termini, remove_nonpolar
        from .receptor_clean import open_prep_model
        res, err = payload
        if err:
            self._fail("✘ " + err, self.h_result)
            return
        try:
            n_seq = pq.copy_seqres(self._h_in_pdb, res["pdb"])      # lets ChimeraX recognise real and cut ends
            if not n_seq:
                self.log_line("Note: no SEQRES records were available, so ChimeraX cannot tell cut chain ends "
                              "apart from real ones; all ends follow PROPKA.")
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            model = open_prep_model(self.session, res["pdb"], self._prep_name)
            removed_c = remove_nonpolar(model)       # polar only; non-polar ones are an option of their own
            self._prep_model = model
            notes = ["PDB2PQR/PROPKA at pH %g" % res["ph"], "%d carbon hydrogens set aside" % removed_c]
            end_notes = fix_termini(self.session, model, self._h_ends, res["table"], res["ph"], self.log_line)
            for r in res["table"]:
                if r["group"] in ("N+", "C-"):
                    n = end_notes.get((r["chain"], r["number"], "N" if r["group"] == "N+" else "C"))
                    if n:
                        r["renamed"] = n
            if res["renamed"]:
                notes.append("renamed: " + ", ".join(res["renamed"][:8]) + (" …" if len(res["renamed"]) > 8 else ""))
            self.pka_table.set_rows(res["table"], res["ph"], "pKa values (saved to %s)" % os.path.basename(res["csv_path"]))
            for r in res["table"]:
                self.log_line("%-16s pKa %6.2f  (model %5.2f)  %+5.2f  %s%s%s" % (
                    r["residue"], r["pka"], r["model_pka"], r["delta"], r["predicted"],
                    ("  [%s]" % r["note"]) if r.get("note") else "",
                    "   <-- within ±%d of the pH" % (1 if r["flag"] == "close" else 2) if r["flag"] else ""))
            flips = res.get("flips") or []
            if flips:
                self.log_line("Flipped by PDB2PQR (hydrogen-bond optimisation): " + ", ".join(
                    "%s (%.1f Å)" % f for f in flips) + ". ChimeraX's own mode never flips. Check them in the scene.")
                notes.append("flipped: " + ", ".join(f[0] for f in flips))
            else:
                self.log_line("PDB2PQR flipped no Asn / Gln / His side chain.")
            self.log_line("Log: %s" % res["log_path"])
            self.log_line("Table: %s" % res["csv_path"])
            self._finish_hydrogens(notes, "pdb2pqr")
        except Exception as err2:
            self._fail("Could not use the PDB2PQR result: %s" % err2, self.h_result)

    # ---- undo / redo of the hydrogen step
    def _swap_model(self, path, which):
        """Replace the copy in the scene by the structure in *path* (idempotent)."""
        from .receptor_clean import open_prep_model
        if getattr(self, "_swap_state", None) == which and self._prep_alive():
            return
        try:
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            self._prep_model = open_prep_model(self.session, path, self._prep_name)
        except Exception as err:
            self._fail("Could not restore %s: %s" % (os.path.basename(path), err), self.h_result)
            return
        self._swap_state = which
        proj = self._project()
        if self._rev is not None and proj is not None:
            if which == "before":
                for st in ("hydrogens", "minimized", "pdbqt"):
                    self._rev["stages"].pop(st, None)
            proj.save()
            self.load_project()
        self._set_undo_state(which)
        self.h_result.setText("Hydrogens: %s." % ("undone, back to the structure before the step" if which == "before"
                                                  else "redone"))
        self.h_result.setStyleSheet(WARN_STYLE)

    def _set_undo_state(self, state):
        self._swap_state = state
        self.h_undo_btn.setEnabled(state == "after")
        self.h_redo_btn.setEnabled(state == "before")

    def _undo_hydrogens(self):
        if self._h_pre and self._swap_state == "after":
            self._swap_model(self._h_pre, "before")

    def _redo_hydrogens(self):
        proj = self._project()
        if proj is None or self._rev is None:
            return
        path = proj.stage_file(self._rev, "hydrogens")
        if self._swap_state == "before" and os.path.isfile(path):
            self._swap_model(path, "after")
            self._rev["stages"]["hydrogens"] = {"file": proj.rel(path)}
            proj.save()
            self.load_project()

    # ------------------------------------------------------- 3 minimize
    def _minimize(self):
        from .receptor_clean import minimize_model
        if self._min_watch is not None:          # running: this button is "Stop"
            self._stop_minimize()
            return
        if not self._prep_alive() and not self._clean_receptor():
            return
        if not (self.ions_cb.isChecked() and self.lig_cb.isChecked()):
            self.log_line("Note: ligands/ions are kept, so minimization may fail for lack of parameters.")
        try:
            self._min_watch = minimize_model(self.session, self._prep_model, self.min_steps.value(),
                                             self.min_conv.isChecked(), on_done=self._minimize_done,
                                             on_energy=self._minimize_energy, log=self.log_line)
        except Exception as err:
            self._fail("Minimization failed: %s" % err)
            return
        self.min_btn.setText("3. Stop minimization")
        self.rec_result.setText("Minimizing ... press Stop to end it early.")
        self.rec_result.setStyleSheet(WARN_STYLE)
        w = self._min_watch
        QTimer.singleShot(8000, lambda: self._check_following(w))

    def _check_following(self, watcher):
        if watcher is self._min_watch and watcher is not None and not watcher.started and not watcher.finished:
            self.log_line("Could not follow the minimization from the status line. Watch ChimeraX's status bar; "
                          "when it says 'Minimization complete', press Stop here and continue.")

    def _minimize_energy(self, step, energy):
        self.log_line("minimization step %d: energy %.2f" % (step, energy))

    def _stop_minimize(self):
        from chimerax.core.commands import run
        w, self._min_watch = self._min_watch, None
        if w is not None:
            w.stop()
        try:
            run(self.session, "minimize stop")
        except Exception:
            pass
        self.min_btn.setText(MIN_TEXT)
        self.rec_result.setText("Minimization stopped. You can run step 4 now.")
        self.rec_result.setStyleSheet(WARN_STYLE)

    def _minimize_done(self, message, energy):
        w, self._min_watch = self._min_watch, None
        self.min_btn.setText(MIN_TEXT)
        self.log_line(message)
        if w is not None and not w.energies:
            self.log_line("(ChimeraX reports energies every 100 steps; runs shorter than that report none.)")
        try:
            self._checkpoint("minimized", energy=energy)
        except Exception as err:
            self.log_line("Could not save the minimized structure: %s" % err)
        self.rec_result.setText("✔ %s. Step 4 converts the minimized model." % message)
        self.rec_result.setStyleSheet(OK_STYLE)
        self.conv_btn.setEnabled(self.controller.prep_ok and not self.controller.prep_busy)

    # ------------------------------------------------------- 4 convert
    def _convert_receptor(self):
        """Step 4: save the cleaned copy and run Meeko."""
        from .receptor_clean import save_model_pdb
        if self._min_watch is not None:
            self.rec_result.setText("Minimization is still running. Wait for it to finish (or press Stop).")
            self.rec_result.setStyleSheet(WARN_STYLE)
            return
        wd = self._workdir()
        if wd is None:
            return
        if not self._prep_alive() and not self._clean_receptor():
            return
        proj = self._project()
        out_dir = proj.rev_dir(self._rev)
        clean_pdb = os.path.join(out_dir, "04_for_meeko.pdb")
        try:
            save_model_pdb(self.session, self._prep_model, clean_pdb)
        except Exception as err:
            self._fail("Saving failed: %s" % err)
            return
        self.rec_result.setText("Running Meeko ...")
        self.rec_result.setStyleSheet("")
        self._convert_rev = self._rev
        self.controller.prepare_receptor_async(clean_pdb, os.path.join(out_dir, self._prep_name), self.bad_cb.isChecked())

    def receptor_done(self, payload):
        res, err = payload
        if err:
            self._fail("✘ " + err)
            return
        path, skipped, rescued = res
        proj, rev = self._project(), getattr(self, "_convert_rev", None)
        if proj is not None and rev is not None:
            proj.set_stage(rev, "pdbqt", path, skipped=list(skipped), rescued=list(rescued))
            proj.data["receptor"]["active"] = rev["id"]
            proj.save()
            self.load_project()
        text = "✔ %s" % path
        if rescued:
            text += "\n✔ Kept after retrying without their hydrogens (Meeko assigned standard hydrogens): %s." % ", ".join(rescued)
        if skipped:
            text += ("\n⚠ Skipped %d residue(s): %s%s. Check that none lie in the binding site (see the log)."
                     % (len(skipped), ", ".join(skipped[:8]), " …" if len(skipped) > 8 else ""))
        self.rec_result.setText(text)
        self.rec_result.setStyleSheet(WARN_STYLE if skipped else OK_STYLE)
        if self._prep_alive():
            try:
                self._prep_model.display = False
            except Exception:
                pass
        text_extra = "\nCheck the result, then press \"Approve for docking\" to offer it on the Inputs tab."
        self.rec_result.setText(self.rec_result.text() + text_extra)
        self._update_approve()
