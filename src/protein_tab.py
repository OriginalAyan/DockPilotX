"""Protein tab: clean a PDB/mmCIF receptor with Dock Prep, set hydrogens, minimize, and convert it with Meeko."""

import os

from Qt.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
    QProgressBar, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)

from .pdbqt_tools import safe_name
from .prep_widgets import BAD_STYLE, OK_STYLE, WARN_STYLE, EnvBanner, PkaTable
from .receptor_clean import SIDE_CHAIN_CHOICES, STANDARDIZE_NOTE

MIN_TEXT = "3. Minimize (optional)"
POLAR_MODES = [
    ("Keep as in the file (add none)", "keep"),
    ("ChimeraX (residue names, else H-bond method)", "chimerax"),
    ("PDB2PQR + PROPKA at pH", "pdb2pqr"),
]


class ProteinTab(QWidget):
    def __init__(self, session, controller, host):
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host          # the tool: workdir_edit, receptor_edit, _receptor_chosen, show_tab
        self._model_ids = []
        self._prep_model = None       # the cleaned copy in the scene
        self._prep_name = None
        self._prep_source = None
        self._min_watch = None
        self._build_ui()
        self.refresh_env()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        lay = QVBoxLayout(self)
        self.banner = EnvBanner(self.controller, self.host)
        lay.addWidget(self.banner)

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
        row.addWidget(QLabel("or open model"))
        self.model_combo = QComboBox()
        self.model_combo.setToolTip("Prepare a copy of a model that is already open. The model itself is not changed.")
        row.addWidget(self.model_combo, 1)
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
        self.rec_btn.setToolTip("Runs Dock Prep on a copy and keeps it in the scene as <name>_prep. Adds no hydrogens.")
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
            "PDB2PQR: PROPKA predicts pKa values and the states follow the pH you enter. Existing hydrogens are replaced.")
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
                                    "Only hydrogens on carbon are added, so polar hydrogens you removed stay removed.")
        gl.addWidget(self.nonpolar_cb)
        self.h_btn = QPushButton("2. Apply hydrogens")
        self.h_btn.setToolTip("Applies the choices above to the cleaned copy (it is cleaned first if needed).")
        self.h_btn.clicked.connect(self._apply_hydrogens)
        gl.addWidget(self.h_btn)
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
                                "Energies go to the ChimeraX Log; this tool tells you when it is done.")
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
        self.bad_cb = QCheckBox("Skip residues Meeko cannot handle")
        self.bad_cb.setToolTip("Residues without a matching template are dropped, with a warning, instead of stopping.")
        self.bad_cb.setChecked(True)
        gl.addWidget(self.bad_cb)
        self.conv_btn = QPushButton("4. Convert with Meeko")
        self.conv_btn.setToolTip("Saves the current <name>_prep model and converts it to PDBQT. "
                                 "If step 1 was not run, it is run first.")
        self.conv_btn.clicked.connect(self._convert_receptor)
        gl.addWidget(self.conv_btn)
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
        """List open atomic models for 'or open model'."""
        self.model_combo.clear()
        self.model_combo.addItem("(use the file above)")
        self._model_ids = [None]
        try:
            from chimerax.atomic import AtomicStructure
            for m in self.session.models.list(type=AtomicStructure):
                self.model_combo.addItem("#%s  %s" % (".".join(str(i) for i in m.id), m.name))
                self._model_ids.append(m)
        except Exception:
            pass

    def _polar_changed(self, *_):
        self.ph_spin.setEnabled(self.polar_combo.currentData() == "pdb2pqr")

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

    def _source(self, out_dir):
        """(source file, name) from the file field or the selected open model; None if neither is usable."""
        idx = self.model_combo.currentIndex()
        if 0 < idx < len(self._model_ids):
            from chimerax.core.commands import run
            m = self._model_ids[idx]
            name = safe_name(m.name) or "receptor"
            source = os.path.join(out_dir, name + "_source.pdb")
            spec = "#" + ".".join(str(i) for i in m.id)
            run(self.session, 'save "%s" %s relModel %s' % (source.replace("\\", "/"), spec, spec))
            return source, name
        source = self.rec_edit.text().strip().strip('"')
        if not source or not os.path.isfile(source):
            return None
        return source, safe_name(os.path.splitext(os.path.basename(source))[0]) or "receptor"

    def _fail(self, text, label=None):
        (label or self.rec_result).setText(text)
        (label or self.rec_result).setStyleSheet(BAD_STYLE)

    def _prep_alive(self):
        m = self._prep_model
        return m is not None and not getattr(m, "deleted", False)

    # ------------------------------------------------------- 1 clean
    def _clean_receptor(self):
        """Step 1: Dock Prep on a copy; the copy stays in the scene."""
        from .receptor_clean import clean_receptor
        wd = self._workdir()
        if wd is None:
            return False
        out_dir = os.path.join(wd, "prepared", "receptor")
        os.makedirs(out_dir, exist_ok=True)
        try:
            src = self._source(out_dir)
            if src is None:
                self._fail("Choose a protein file or an open model.")
                return False
            source, name = src
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            model, summary = clean_receptor(
                self.session, source, name, chains=self.chains_edit.text(),
                del_water=self.water_cb.isChecked(), del_ions=self.ions_cb.isChecked(),
                del_ligands=self.lig_cb.isChecked(), del_altlocs=self.alt_cb.isChecked(),
                del_incomplete=self.bb_cb.isChecked(), standardize=self.std_cb.isChecked(),
                side_chains=self.sc_combo.currentData() if self.sc_cb.isChecked() else "", log=self.log_line)
        except Exception as err:
            self._fail("Cleaning failed: %s" % err)
            return False
        self._prep_model, self._prep_name, self._prep_source = model, name, source
        self.pka_table.setVisible(False)
        self.log_line(summary)
        self.rec_result.setText("✔ %s Step 2 sets hydrogens (optional)." % summary)
        self.rec_result.setStyleSheet(OK_STYLE)
        return True

    # ------------------------------------------------------- 2 hydrogens
    def _apply_hydrogens(self):
        from .hydrogens import add_hydrogens, count_hydrogens
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
        if mode == "pdb2pqr":
            from .receptor_clean import save_model_pdb
            out_dir = os.path.join(wd, "prepared", "receptor", "pdb2pqr")
            os.makedirs(out_dir, exist_ok=True)
            in_pdb = os.path.join(out_dir, self._prep_name + "_before_pdb2pqr.pdb")
            try:
                save_model_pdb(self.session, self._prep_model, in_pdb)
            except Exception as err:
                self._fail("Saving failed: %s" % err, self.h_result)
                return
            self.h_result.setText("Running PDB2PQR + PROPKA ...")
            self.h_result.setStyleSheet("")
            self.controller.run_pdb2pqr_async(in_pdb, out_dir, self._prep_name, self.ph_spin.value())
            return
        notes = []
        try:
            if mode == "chimerax":
                added, wrong = add_hydrogens(self.session, self._prep_model, "polar", self.log_line)
                notes.append("%d polar hydrogens added" % added)
            self._finish_hydrogens(notes)
        except Exception as err:
            self._fail("Adding hydrogens failed: %s" % err, self.h_result)

    def _finish_hydrogens(self, notes):
        """Non-polar step (if wanted) and the summary line."""
        from .hydrogens import add_hydrogens, count_hydrogens
        if self.nonpolar_cb.isChecked():
            added, wrong = add_hydrogens(self.session, self._prep_model, "nonpolar", self.log_line)
            notes.append("%d non-polar hydrogens added" % added)
        pol, nonpol = count_hydrogens(self._prep_model)
        notes.append("now %d polar and %d non-polar hydrogens" % (pol, nonpol))
        text = "✔ " + "; ".join(notes) + "."
        self.log_line(text)
        self.h_result.setText(text)
        self.h_result.setStyleSheet(OK_STYLE)

    def pdb2pqr_done(self, payload):
        from .hydrogens import remove_nonpolar
        from .receptor_clean import open_prep_model
        res, err = payload
        if err:
            self._fail("✘ " + err, self.h_result)
            return
        try:
            if self._prep_alive():
                self.session.models.close([self._prep_model])
            model = open_prep_model(self.session, res["pdb"], self._prep_name)
            n_c = remove_nonpolar(model)       # polar only; non-polar ones are an option of their own
            self._prep_model = model
            notes = ["PDB2PQR/PROPKA at pH %g" % res["ph"]]
            if res["renamed"]:
                notes.append("renamed: " + ", ".join(res["renamed"][:8]) + (" …" if len(res["renamed"]) > 8 else ""))
            self.pka_table.set_rows(res["table"], res["ph"], "pKa values (saved to %s)" % os.path.basename(res["csv_path"]))
            for r in res["table"]:
                self.log_line("%-16s pKa %6.2f  (model %5.2f)  %+5.2f  %s%s" % (
                    r["residue"], r["pka"], r["model_pka"], r["delta"], r["predicted"],
                    "   <-- within ±%d of the pH" % (1 if r["flag"] == "close" else 2) if r["flag"] else ""))
            self.log_line("Log: %s" % res["log_path"])
            self.log_line("Table: %s" % res["csv_path"])
            self._finish_hydrogens(notes)
        except Exception as err2:
            self._fail("Could not use the PDB2PQR result: %s" % err2, self.h_result)

    # ------------------------------------------------------- 3 minimize
    def _minimize(self):
        from .hydrogens import lacks_nonpolar
        from .receptor_clean import minimize_model
        if self._min_watch is not None:          # running: this button is "Stop"
            self._stop_minimize()
            return
        if not self._prep_alive() and not self._clean_receptor():
            return
        if not (self.ions_cb.isChecked() and self.lig_cb.isChecked()):
            self.log_line("Note: ligands/ions are kept, so minimization may fail for lack of parameters.")
        if lacks_nonpolar(self._prep_model):
            self.log_line("Note: the model has no non-polar hydrogens. Minimization works best with all hydrogens "
                          "(step 2: tick 'Add non-polar hydrogens').")
        try:
            self._min_watch = minimize_model(self.session, self._prep_model, self.min_steps.value(),
                                             self.min_conv.isChecked(), on_done=self._minimize_done)
        except Exception as err:
            self._fail("Minimization failed: %s" % err)
            return
        self.min_btn.setText("3. Stop minimization")
        self.rec_result.setText("Minimizing ... energies are written to the ChimeraX Log.")
        self.rec_result.setStyleSheet(WARN_STYLE)

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
        self._min_watch = None
        self.min_btn.setText(MIN_TEXT)
        self.rec_result.setText("✔ %s. Step 4 converts the minimized model." % message)
        self.rec_result.setStyleSheet(OK_STYLE)
        self.log_line(message)
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
        out_dir = os.path.join(wd, "prepared", "receptor")
        clean_pdb = os.path.join(out_dir, self._prep_name + "_clean.pdb")
        try:
            save_model_pdb(self.session, self._prep_model, clean_pdb)
        except Exception as err:
            self._fail("Saving failed: %s" % err)
            return
        self.rec_result.setText("Running Meeko ...")
        self.rec_result.setStyleSheet("")
        self.controller.prepare_receptor_async(clean_pdb, os.path.join(out_dir, self._prep_name), self.bad_cb.isChecked())

    def receptor_done(self, payload):
        res, err = payload
        if err:
            self._fail("✘ " + err)
            return
        path, skipped, rescued = res
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
        self.host.receptor_edit.setText(path)
        self.host._receptor_chosen(open_model=True)
