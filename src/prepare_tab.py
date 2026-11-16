"""Prepare tab: clean a PDB/mmCIF receptor in ChimeraX and convert receptor and ligands to PDBQT with Meeko."""

import os
import re
import tempfile

from Qt.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QListWidget, QPlainTextEdit,
    QProgressBar, QPushButton, QVBoxLayout, QWidget,
)

from . import prep
from .pdbqt_tools import safe_name

OK_STYLE = "color: #2e7d32; font-weight: bold;"
BAD_STYLE = "color: #c62828; font-weight: bold;"
WARN_STYLE = "color: #b26a00;"


class PrepareTab(QWidget):
    def __init__(self, session, controller, host):
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host          # the tool: workdir_edit, receptor_edit, _receptor_chosen, _add_ligand_paths, show_tab
        self._model_ids = []
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
        self.water_cb = QCheckBox("Remove water")
        self.ions_cb = QCheckBox("Remove ions / metals")
        self.lig_cb = QCheckBox("Remove ligands / cofactors")
        self.addh_cb = QCheckBox("Add hydrogens (ChimeraX)")
        self.bad_cb = QCheckBox("Skip residues Meeko cannot handle")
        self.addh_cb.setToolTip("Meeko needs hydrogens on the receptor.")
        self.bad_cb.setToolTip("Residues without a matching template are dropped, with a warning, instead of stopping.")
        for cb in (self.water_cb, self.ions_cb, self.lig_cb, self.addh_cb, self.bad_cb):
            cb.setChecked(True)
            gl.addWidget(cb)
        self.rec_btn = QPushButton("Prepare receptor")
        self.rec_btn.clicked.connect(self._prepare_receptor)
        gl.addWidget(self.rec_btn)
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
        for b in (self.rec_btn, self.lig_btn):
            b.setEnabled(ok)
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
            self.rec_result.setText("Choose a working directory on the Settings tab first.")
            self.rec_result.setStyleSheet(BAD_STYLE)
            self.lig_result.setText("Choose a working directory on the Settings tab first.")
            self.lig_result.setStyleSheet(BAD_STYLE)
            return None
        return os.path.abspath(w)

    def _browse_receptor(self):
        p, _ = QFileDialog.getOpenFileName(self, "Select a protein structure", self.rec_edit.text(),
                                           "Structures (*.pdb *.ent *.cif *.mmcif);;All files (*)")
        if p:
            self.rec_edit.setText(p)
            self.model_combo.setCurrentIndex(0)

    def _prepare_receptor(self):
        from .receptor_clean import clean_receptor
        wd = self._workdir()
        if wd is None:
            return
        out_dir = os.path.join(wd, "prepared", "receptor")
        os.makedirs(out_dir, exist_ok=True)
        idx = self.model_combo.currentIndex()
        source, name = self.rec_edit.text().strip().strip('"'), None
        try:
            if idx > 0 and idx < len(self._model_ids):
                from chimerax.core.commands import run
                m = self._model_ids[idx]
                name = safe_name(m.name) or "receptor"
                source = os.path.join(out_dir, name + "_source.pdb")
                spec = "#" + ".".join(str(i) for i in m.id)
                run(self.session, 'save "%s" %s relModel %s' % (source.replace("\\", "/"), spec, spec))
            elif not source or not os.path.isfile(source):
                self.rec_result.setText("Choose a protein file or an open model.")
                self.rec_result.setStyleSheet(BAD_STYLE)
                return
            else:
                name = safe_name(os.path.splitext(os.path.basename(source))[0]) or "receptor"
            clean_pdb = os.path.join(out_dir, name + "_clean.pdb")
            summary = clean_receptor(self.session, source, clean_pdb, chains=self.chains_edit.text(),
                                     remove_water=self.water_cb.isChecked(), remove_ions=self.ions_cb.isChecked(),
                                     remove_ligands=self.lig_cb.isChecked(), add_h=self.addh_cb.isChecked(),
                                     log=self.log_line)
        except Exception as err:
            self.rec_result.setText("Cleaning failed: %s" % err)
            self.rec_result.setStyleSheet(BAD_STYLE)
            return
        self.log_line(summary)
        self.rec_result.setText(summary + " Running Meeko ...")
        self.rec_result.setStyleSheet("")
        self.controller.prepare_receptor_async(clean_pdb, os.path.join(out_dir, name), self.bad_cb.isChecked())

    def receptor_done(self, payload):
        res, err = payload
        if err:
            self.rec_result.setText("✘ " + err)
            self.rec_result.setStyleSheet(BAD_STYLE)
            return
        path, warnings = res
        text = "✔ %s" % path
        if warnings:
            text += "\n⚠ Skipped %d residue(s): %s%s (see the log)." % (len(warnings), ", ".join(warnings[:8]), " …" if len(warnings) > 8 else "")
        self.rec_result.setText(text)
        self.rec_result.setStyleSheet(WARN_STYLE if warnings else OK_STYLE)
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
