"""Prepare tab: clean a PDB/mmCIF receptor in ChimeraX and convert receptor and ligands to PDBQT with Meeko."""

import os
import re
import tempfile

from Qt.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget, QPlainTextEdit,
    QProgressBar, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)

from . import prep
from .receptor_clean import SIDE_CHAIN_CHOICES, STANDARDIZE_NOTE
from .pdbqt_tools import safe_name

OK_STYLE = "color: #2e7d32; font-weight: bold;"
BAD_STYLE = "color: #c62828; font-weight: bold;"
WARN_STYLE = "color: #b26a00;"
MIN_TEXT = "2. Minimize (optional)"


class PrepareTab(QWidget):
    def __init__(self, session, controller, host):
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host          # the tool: workdir_edit, receptor_edit, _receptor_chosen, _add_ligand_paths, show_tab
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
        self.env_label = QLabel("")
        self.env_label.setWordWrap(True)
        lay.addWidget(self.env_label)
        goto = QPushButton("Set up Meeko on the Settings tab")
        goto.clicked.connect(lambda: self.host.show_tab(0))
        self.goto_btn = goto
        lay.addWidget(goto)

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
                note = QLabel("\n".join("   \u2022 " + t for t in STANDARDIZE_NOTE))
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
        self.h_cb = QCheckBox("Add hydrogens")
        self.h_cb.setChecked(True)
        self.h_cb.setToolTip("Dock Prep option. Meeko adds standard polar hydrogens itself, so unchecking is possible,\n"
                             "but then every histidine gets Meeko's default protonation.")
        gl.addWidget(self.h_cb)
        self.h_only_cb = QCheckBox("Only if the structure has none")
        self.h_only_cb.setChecked(True)
        self.h_only_cb.setToolTip("Checked: a structure that already has hydrogens is left alone, so protonation you set\n"
                                  "earlier (e.g. HID / HIE / HIP) is kept.\n"
                                  "Unchecked: Dock Prep also adds hydrogens that are missing from a partly protonated structure.")
        self.h_only_cb.setStyleSheet("margin-left: 18px;")
        self.h_cb.toggled.connect(self.h_only_cb.setEnabled)
        gl.addWidget(self.h_only_cb)
        self.bad_cb = QCheckBox("Skip residues Meeko cannot handle")
        self.bad_cb.setToolTip("Residues without a matching template are dropped, with a warning, instead of stopping.")
        self.bad_cb.setChecked(True)
        gl.addWidget(self.bad_cb)
        self.rec_btn = QPushButton("1. Clean receptor")
        self.rec_btn.setToolTip("Runs Dock Prep on a copy and keeps it in the scene as <name>_prep.")
        self.rec_btn.clicked.connect(self._clean_receptor)
        gl.addWidget(self.rec_btn)
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
        self.conv_btn = QPushButton("3. Convert with Meeko")
        self.conv_btn.setToolTip("Saves the current <name>_prep model and converts it to PDBQT. "
                                 "If step 1 was not run, it is run first.")
        self.conv_btn.clicked.connect(self._convert_receptor)
        gl.addWidget(self.conv_btn)
        self.rec_result = QLabel("")
        self.rec_result.setWordWrap(True)
        gl.addWidget(self.rec_result)
        lay.addWidget(g)

        g = QGroupBox("Ligands (.sdf .mol .mol2 .smi)")
        gl = QVBoxLayout(g)
        row = QHBoxLayout()
        for text, slot in (("Add files…", self._add_files), ("Add folder…", self._add_folder),
                           ("Remove selected", self._remove), ("Clear", lambda: self.lig_list.clear())):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        gl.addLayout(row)
        self.lig_list = QListWidget()
        self.lig_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        self.lig_list.setMinimumHeight(90)
        gl.addWidget(self.lig_list)
        self.lig_h_cb = QCheckBox("Add hydrogens and 3D coordinates if missing")
        self.lig_h_cb.setToolTip("Done with RDKit. Existing 3D structures with hydrogens are used as they are.\n"
                                 "Protonation states are not changed.")
        self.lig_h_cb.setChecked(True)
        gl.addWidget(self.lig_h_cb)
        self.lig_add_cb = QCheckBox("Add the results to the Inputs tab")
        self.lig_add_cb.setChecked(True)
        gl.addWidget(self.lig_add_cb)
        self.lig_btn = QPushButton("Prepare ligands")
        self.lig_btn.clicked.connect(self._prepare_ligands)
        gl.addWidget(self.lig_btn)
        self.lig_result = QLabel("")
        self.lig_result.setWordWrap(True)
        gl.addWidget(self.lig_result)
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
        c = self.controller
        if c.prep_ok:
            self.env_label.setText("✔ " + (c.prep_message or "Meeko environment ready."))
            self.env_label.setStyleSheet(OK_STYLE)
        else:
            self.env_label.setText("✘ " + (c.prep_message or "Meeko is not set up yet."))
            self.env_label.setStyleSheet(BAD_STYLE)
        self.goto_btn.setVisible(not c.prep_ok)
        self.set_busy(c.prep_busy)

    def set_busy(self, busy):
        ok = self.controller.prep_ok and not busy
        for b in (self.conv_btn, self.lig_btn):
            b.setEnabled(ok)
        self.rec_btn.setEnabled(not busy)
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

    # ------------------------------------------------------------ receptor
    def _workdir(self):
        w = self.host.workdir_edit.text().strip().strip('"')
        if not w:
            for lab in (self.rec_result, self.lig_result):
                lab.setText("Choose a working directory on the Settings tab first.")
                lab.setStyleSheet(BAD_STYLE)
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

    def _fail(self, text):
        self.rec_result.setText(text)
        self.rec_result.setStyleSheet(BAD_STYLE)

    def _prep_alive(self):
        m = self._prep_model
        return m is not None and not getattr(m, "deleted", False)

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
                side_chains=self.sc_combo.currentData() if self.sc_cb.isChecked() else "",
                hydrogens=self._hydrogen_mode(), log=self.log_line)
        except Exception as err:
            self._fail("Cleaning failed: %s" % err)
            return False
        self._prep_model, self._prep_name, self._prep_source = model, name, source
        self.log_line(summary)
        self.rec_result.setText("✔ %s Step 2 (minimize) is optional; step 3 converts it." % summary)
        self.rec_result.setStyleSheet(OK_STYLE)
        return True

    def _hydrogen_mode(self):
        if not self.h_cb.isChecked():
            return "never"
        return "auto" if self.h_only_cb.isChecked() else "always"

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
                                             self.min_conv.isChecked(), on_done=self._minimize_done)
        except Exception as err:
            self._fail("Minimization failed: %s" % err)
            return
        self.min_btn.setText("2. Stop minimization")
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
        self.rec_result.setText("Minimization stopped. You can run step 3 now.")
        self.rec_result.setStyleSheet(WARN_STYLE)

    def _minimize_done(self, message, energy):
        self._min_watch = None
        self.min_btn.setText(MIN_TEXT)
        self.rec_result.setText("\u2714 %s. Step 3 converts the minimized model." % message)
        self.rec_result.setStyleSheet(OK_STYLE)
        self.log_line(message)
        self.conv_btn.setEnabled(self.controller.prep_ok and not self.controller.prep_busy)

    def _convert_receptor(self):
        """Step 3: save the cleaned copy and run Meeko."""
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

    # ------------------------------------------------------------- ligands
    def _add_files(self):
        ps, _ = QFileDialog.getOpenFileNames(self, "Select ligand files", "",
                                             "Ligands (*.sdf *.mol *.mol2 *.smi *.smiles);;All files (*)")
        self._add(ps)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Select a folder with ligand files")
        if d:
            self._add([d])

    def _add(self, paths):
        have = {self.lig_list.item(i).text() for i in range(self.lig_list.count())}
        for p in paths:
            if p and p not in have:
                self.lig_list.addItem(p)

    def _remove(self):
        for it in self.lig_list.selectedItems():
            self.lig_list.takeItem(self.lig_list.row(it))

    def _prepare_ligands(self):
        wd = self._workdir()
        if wd is None:
            return
        inputs = [self.lig_list.item(i).text() for i in range(self.lig_list.count())]
        if not prep.list_ligand_files(inputs):
            self.lig_result.setText("Add ligand files first (.sdf, .mol, .mol2, .smi).")
            self.lig_result.setStyleSheet(BAD_STYLE)
            return
        base = os.path.join(wd, "prepared")
        self.lig_result.setText("Working ...")
        self.lig_result.setStyleSheet("")
        self.bar.setRange(0, 100)
        self.controller.prepare_ligands_async(inputs, os.path.join(base, "ligands"), os.path.join(base, "staging"),
                                              self.lig_h_cb.isChecked())

    def progress(self, payload):
        done, total = payload
        if total:
            self.bar.setRange(0, 100)
            self.bar.setValue(int(100 * done / total))

    def ligands_done(self, payload):
        res, err = payload
        if err:
            self.lig_result.setText("✘ " + err)
            self.lig_result.setStyleSheet(BAD_STYLE)
            return
        n, bad = len(res["files"]), res["failed"]
        text = "✔ %d ligand(s) prepared in prepared/ligands." % n if n else "✘ No ligands were prepared."
        if bad:
            text += "\n⚠ %d failed: %s%s" % (len(bad), ", ".join(b[0] for b in bad[:5]), " …" if len(bad) > 5 else "")
        self.lig_result.setText(text)
        self.lig_result.setStyleSheet(OK_STYLE if n and not bad else (WARN_STYLE if n else BAD_STYLE))
        if n and self.lig_add_cb.isChecked():
            self.host._add_ligand_paths(res["files"])
