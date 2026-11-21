"""Ligands tab: convert ligand files (.sdf .mol .mol2 .smi) to PDBQT with RDKit staging and Meeko."""

import os

from Qt.QtWidgets import (
    QCheckBox, QFileDialog, QGroupBox, QHBoxLayout, QLabel, QListWidget, QPlainTextEdit, QProgressBar, QPushButton,
    QVBoxLayout, QWidget,
)

from . import prep
from .prep_widgets import BAD_STYLE, OK_STYLE, WARN_STYLE, EnvBanner


class LigandTab(QWidget):
    def __init__(self, session, controller, host):
        super().__init__()
        self.session = session
        self.controller = controller
        self.host = host
        self._build_ui()
        self.refresh_env()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        self.banner = EnvBanner(self.controller, self.host)
        lay.addWidget(self.banner)

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

    def refresh_env(self):
        self.banner.refresh()
        self.set_busy(self.controller.prep_busy)

    def set_busy(self, busy):
        self.lig_btn.setEnabled(self.controller.prep_ok and not busy)
        self.cancel_btn.setVisible(busy)
        self.bar.setVisible(busy)
        if busy:
            self.bar.setRange(0, 0)

    def log_line(self, text):
        self.log.appendPlainText(str(text))

    def _workdir(self):
        w = self.host.workdir_edit.text().strip().strip('"')
        if not w:
            self.lig_result.setText("Choose a working directory on the Settings tab first.")
            self.lig_result.setStyleSheet(BAD_STYLE)
            return None
        return os.path.abspath(w)

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
