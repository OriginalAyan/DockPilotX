"""Small widgets shared by the Protein and Ligands tabs."""

from Qt.QtCore import Qt
from Qt.QtGui import QBrush, QColor
from Qt.QtWidgets import (
    QAbstractItemView, QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

OK_STYLE = "color: #2e7d32; font-weight: bold;"
BAD_STYLE = "color: #c62828; font-weight: bold;"
WARN_STYLE = "color: #b26a00;"


def is_dark(widget):
    """True when the widget's palette is a dark theme (ChimeraX has light and dark modes)."""
    try:
        return widget.palette().window().color().lightness() < 128
    except Exception:
        return False


def flag_colors(dark):
    """(strong, soft, text) colours for the two pKa flags. Hue tells them apart (red-orange vs yellow), so they
    stay distinct in both themes; a dark theme uses dark fills with white text."""
    if dark:
        return QColor("#9b3328"), QColor("#6b5b16"), QColor("#ffffff")
    return QColor("#f4a39a"), QColor("#fbe9a6"), QColor("#000000")


class EnvBanner(QWidget):
    """One line saying whether Meeko is ready, with a shortcut to the Settings tab."""

    def __init__(self, controller, host, what="Meeko"):
        super().__init__()
        self.controller, self.host, self.what = controller, host, what
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel("")
        self.label.setWordWrap(True)
        lay.addWidget(self.label)
        self.btn = QPushButton("Set up %s on the Settings tab" % what)
        self.btn.clicked.connect(lambda: self.host.show_tab(0))
        lay.addWidget(self.btn)

    def refresh(self):
        c = self.controller
        if c.prep_ok:
            self.label.setText("✔ " + (c.prep_message or "Meeko environment ready."))
            self.label.setStyleSheet(OK_STYLE)
        else:
            self.label.setText("✘ " + (c.prep_message or "Meeko is not set up yet."))
            self.label.setStyleSheet(BAD_STYLE)
        self.btn.setVisible(not c.prep_ok)


class PkaTable(QWidget):
    """PROPKA results: rows within 1 pH unit of the pH get a strong colour, within 2 a soft one."""

    HEADERS = ["Residue", "pKa", "Model pKa", "pKa − pH", "Predicted", "Applied", "Check"]

    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.title = QLabel("")
        self.title.setWordWrap(True)
        lay.addWidget(self.title)
        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setMinimumHeight(150)
        lay.addWidget(self.table)
        row = QHBoxLayout()
        self.strong_lbl, self.soft_lbl = QLabel(" within ±1 of the pH "), QLabel(" within ±2 of the pH ")
        row.addWidget(self.strong_lbl)
        row.addWidget(self.soft_lbl)
        row.addStretch(1)
        lay.addLayout(row)
        self.note = QLabel("Manual verification may be needed where the pKa is close to the pH. "
                           "The states were applied as PDB2PQR returned them.")
        self.note.setWordWrap(True)
        lay.addWidget(self.note)
        self._rows = []
        self.setVisible(False)
        self._paint_legend()

    def _paint_legend(self):
        strong, soft, text = flag_colors(is_dark(self))
        for lbl, col in ((self.strong_lbl, strong), (self.soft_lbl, soft)):
            lbl.setStyleSheet("background-color: %s; color: %s; padding: 1px 4px;" % (col.name(), text.name()))

    def set_rows(self, rows, ph, title=""):
        self._rows = list(rows)
        self._ph = ph
        self.title.setText(title)
        self.table.setRowCount(0)
        for r in self._rows:
            i = self.table.rowCount()
            self.table.insertRow(i)
            vals = [r["residue"], "%.2f" % r["pka"], "%.2f" % r["model_pka"], "%+.2f" % r["delta"], r["predicted"],
                    r["renamed"] or "-", {"close": "±1", "near": "±2"}.get(r["flag"], "")]
            for j, v in enumerate(vals):
                it = QTableWidgetItem(v)
                it.setTextAlignment(Qt.AlignmentFlag.AlignCenter if j else Qt.AlignmentFlag.AlignLeft)
                self.table.setItem(i, j, it)
        self.table.resizeColumnsToContents()
        self.recolor()
        self.setVisible(bool(self._rows))

    def recolor(self):
        strong, soft, text = flag_colors(is_dark(self))
        self._paint_legend()
        for i, r in enumerate(self._rows):
            col = strong if r["flag"] == "close" else (soft if r["flag"] == "near" else None)
            for j in range(self.table.columnCount()):
                it = self.table.item(i, j)
                if it is None:
                    continue
                if col is None:
                    it.setData(Qt.ItemDataRole.BackgroundRole, None)
                    it.setData(Qt.ItemDataRole.ForegroundRole, None)
                else:
                    it.setBackground(QBrush(col))
                    it.setForeground(QBrush(text))

    def changeEvent(self, event):
        super().changeEvent(event)
        try:
            from Qt.QtCore import QEvent
            if event.type() in (QEvent.Type.PaletteChange, QEvent.Type.StyleChange):
                self.recolor()
        except Exception:
            pass

    def showEvent(self, event):
        super().showEvent(event)
        self.recolor()      # the theme may have changed while the tab was hidden
