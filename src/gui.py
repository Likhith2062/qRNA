import csv
import sys
import tempfile
from pathlib import Path

import main
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QFrame, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QMainWindow, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QVBoxLayout, QWidget, QComboBox
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def draw_structure(sequence, dot_bracket, filename, title):
    """Build a simple RNA base-pair arc diagram from dot-bracket notation."""
    n = len(sequence)
    fig, ax = plt.subplots(figsize=(max(8, min(18, n * 0.38)), 4.8))
    ax.plot(range(n), [0] * n, linewidth=1.4)

    for i, base in enumerate(sequence):
        ax.text(i, -0.08, base, ha="center", va="top",
                fontsize=10, fontweight="bold")

    stack, pairs = [], []
    for i, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            pairs.append((stack.pop(), i))

    import numpy as np
    for left, right in pairs:
        center = (left + right) / 2
        radius = (right - left) / 2
        theta = np.linspace(0, np.pi, 80)
        ax.plot(center + radius * np.cos(theta),
                radius * 0.72 * np.sin(theta) + 0.12,
                linewidth=1.5)

    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlim(-1, max(n, 1))
    max_arc = max([((b - a) / 2) * 0.72 for a, b in pairs], default=1)
    ax.set_ylim(-0.45, max(1, max_arc + 0.5))
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


class PredictionWorker(QThread):
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, sequence):
        super().__init__()
        self.sequence = sequence

    def run(self):
        try:
            # The GUI deliberately delegates all scientific computation to main.predict().
            self.finished.emit(main.predict(self.sequence))
        except Exception as exc:
            self.error.emit(str(exc))


class BatchWorker(QThread):
    progress = Signal(int)
    row_finished = Signal(int, object, str)
    finished = Signal()

    def __init__(self, sequences):
        super().__init__()
        self.sequences = sequences

    def run(self):
        for i, sequence in enumerate(self.sequences):
            try:
                result = main.predict(sequence)
                self.row_finished.emit(i, result, "")
            except Exception as exc:
                self.row_finished.emit(i, None, str(exc))
            self.progress.emit(int((i + 1) * 100 / len(self.sequences)))
        self.finished.emit()


class StructureCard(QFrame):
    def __init__(self, title, structure, image):
        super().__init__()
        self.setObjectName("StructureCard")
        layout = QVBoxLayout(self)

        title_label = QLabel(title)
        title_label.setObjectName("CardTitle")

        structure_label = QLabel(structure)
        structure_label.setObjectName("StructureText")
        structure_label.setAlignment(Qt.AlignCenter)
        structure_label.setTextInteractionFlags(Qt.TextSelectableByMouse)

        image_label = QLabel()
        image_label.setAlignment(Qt.AlignCenter)
        pixmap = QPixmap(image)
        image_label.setPixmap(
            pixmap.scaled(850, 300, Qt.KeepAspectRatio,
                          Qt.SmoothTransformation)
        )

        layout.addWidget(title_label)
        layout.addWidget(structure_label)
        layout.addWidget(image_label)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.batch_results = []
        self.temp_dir = Path(tempfile.mkdtemp(prefix="rna_predictor_"))

        self.setWindowTitle("RNA Secondary Structure Predictor")
        self.resize(1280, 850)

        self.tabs = QTabWidget()
        self.tabs.addTab(self.single_tab(), "Single Sequence")
        self.tabs.addTab(self.batch_tab(), "Batch Analysis")
        self.setCentralWidget(self.tabs)
        self.apply_style()

    def single_tab(self):
        page = QWidget()
        root = QVBoxLayout(page)

        h = QLabel("RNA Secondary Structure Predictor")
        h.setObjectName("Header")
        s = QLabel("QUBO-based folding compared with ViennaRNA MFE")
        s.setObjectName("Subtitle")
        root.addWidget(h)
        root.addWidget(s)

        box = QGroupBox("RNA Sequence")
        layout = QVBoxLayout(box)
        self.sequence = QTextEdit()
        self.sequence.setPlaceholderText(
            "Enter RNA sequence, e.g. GGCGAAAUCGCCUUUGGCGAAAUCGCC"
        )
        self.sequence.setMaximumHeight(100)
        self.predict = QPushButton("Predict Structure")
        self.predict.clicked.connect(self.run_single)
        layout.addWidget(self.sequence)
        layout.addWidget(self.predict)
        root.addWidget(box)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.hide()
        root.addWidget(self.progress)

        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.results = QWidget()
        self.results_layout = QVBoxLayout(self.results)
        label = QLabel("Enter an RNA sequence and click Predict Structure.")
        label.setAlignment(Qt.AlignCenter)
        self.results_layout.addWidget(label)
        self.scroll.setWidget(self.results)
        root.addWidget(self.scroll)
        return page

    def batch_tab(self):
        page = QWidget()
        root = QVBoxLayout(page)

        h = QLabel("Batch RNA Analysis")
        h.setObjectName("Header")
        root.addWidget(h)
        root.addWidget(QLabel(
            "Upload a CSV or Excel file containing one RNA sequence per row."
        ))

        row = QHBoxLayout()
        browse = QPushButton("Upload CSV / Excel")
        browse.clicked.connect(self.load_file)
        self.file_label = QLabel("No file selected")
        self.columns = QComboBox()
        self.columns.setEnabled(False)
        self.batch_run = QPushButton("Run Batch Prediction")
        self.batch_run.setEnabled(False)
        self.batch_run.clicked.connect(self.run_batch)

        row.addWidget(browse)
        row.addWidget(self.file_label, 1)
        row.addWidget(QLabel("Sequence column:"))
        row.addWidget(self.columns)
        row.addWidget(self.batch_run)
        root.addLayout(row)

        self.batch_progress = QProgressBar()
        root.addWidget(self.batch_progress)

        self.table = QTableWidget()
        self.table.setColumnCount(11)
        self.table.setHorizontalHeaderLabels([
            "ID", "Sequence", "Stems", "QUBO Coefficients",
            "QUBO Structure", "ViennaRNA", "TP", "TN", "FP", "FN", "MCC"
        ])
        self.table.setAlternatingRowColors(True)
        self.table.cellDoubleClicked.connect(self.batch_detail)
        root.addWidget(self.table)

        export = QPushButton("Export Results as CSV")
        export.clicked.connect(self.export)
        root.addWidget(export)

        self.batch_data = []
        return page

    def run_single(self):
        sequence = "".join(self.sequence.toPlainText().split()).upper()
        if not sequence:
            QMessageBox.warning(self, "Missing sequence", "Please enter an RNA sequence.")
            return
        self.predict.setEnabled(False)
        self.progress.show()
        self.worker = PredictionWorker(sequence)
        self.worker.finished.connect(self.show_result)
        self.worker.error.connect(self.single_error)
        self.worker.start()

    def single_error(self, message):
        self.predict.setEnabled(True)
        self.progress.hide()
        QMessageBox.critical(self, "Prediction Error", message)

    def clear_results(self):
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

    def show_result(self, result):
        self.predict.setEnabled(True)
        self.progress.hide()
        self.clear_results()

        grid = QGridLayout()
        stats = [
            ("Sequence Length", len(result.sequence)),
            ("Candidate Stems", len(result.stems)),
            ("QUBO Coefficients", len(result.qubo)),
            ("ViennaRNA MFE", f"{result.mfe:.2f} kcal/mol"),
            ("QUBO Energy", f"{result.solver_energy:.4f}"),
            ("TP", result.metrics["TP"]),
            ("TN", result.metrics["TN"]),
            ("FP", result.metrics["FP"]),
            ("FN", result.metrics["FN"]),
            ("MCC", f"{result.metrics['MCC']:.4f}"),
        ]
        for i, (name, value) in enumerate(stats):
            card = QFrame()
            card.setObjectName("StatCard")
            l = QVBoxLayout(card)
            n = QLabel(name); n.setObjectName("StatName")
            v = QLabel(str(value)); v.setObjectName("StatValue")
            l.addWidget(n); l.addWidget(v)
            grid.addWidget(card, i // 5, i % 5)
        self.results_layout.addLayout(grid)

        group = QGroupBox("Structure Comparison")
        row = QHBoxLayout(group)
        qimg = self.make_image(result.sequence, result.prediction, "QUBO Prediction")
        vimg = self.make_image(result.sequence, result.reference, "ViennaRNA MFE")
        row.addWidget(StructureCard("QUBO Predicted Structure", result.prediction, qimg))
        row.addWidget(StructureCard("ViennaRNA MFE Structure", result.reference, vimg))
        self.results_layout.addWidget(group)

        note = QLabel(
            "TP/TN/FP/FN are adjacency-matrix entry counts. "
            "MCC summarizes the prediction against the ViennaRNA reference."
        )
        note.setWordWrap(True)
        self.results_layout.addWidget(note)
        self.results_layout.addStretch()

    def make_image(self, sequence, structure, title):
        path = self.temp_dir / (
            title.lower().replace(" ", "_") + f"_{abs(hash(structure))}.png"
        )
        draw_structure(sequence, structure, str(path), title)
        return str(path)

    def load_file(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Select RNA sequence file", "", "Data Files (*.csv *.xlsx *.xls)"
        )
        if not filename:
            return
        try:
            suffix = Path(filename).suffix.lower()
            if suffix == ".csv":
                with open(filename, "r", encoding="utf-8-sig", newline="") as f:
                    reader = csv.DictReader(f)
                    headers = reader.fieldnames or []
                    rows = list(reader)
            else:
                import pandas as pd
                df = pd.read_excel(filename).fillna("")
                headers = [str(c) for c in df.columns]
                rows = df.to_dict("records")
            self.batch_data = rows
            self.file_label.setText(Path(filename).name)
            self.columns.clear()
            self.columns.addItems(headers)
            self.columns.setEnabled(True)
            self.batch_run.setEnabled(True)
        except Exception as exc:
            QMessageBox.critical(self, "File Error", str(exc))

    def run_batch(self):
        col = self.columns.currentText()
        sequences = ["".join(str(r.get(col, "")).split()).upper() for r in self.batch_data]
        sequences = [s for s in sequences if s]
        if not sequences:
            QMessageBox.warning(self, "No sequences", "No non-empty RNA sequences were found.")
            return

        self.batch_results = [None] * len(sequences)
        self.table.setRowCount(len(sequences))
        for i, s in enumerate(sequences):
            self.table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
            self.table.setItem(i, 1, QTableWidgetItem(s))

        self.batch_run.setEnabled(False)
        self.batch_progress.setValue(0)
        self.batch_worker = BatchWorker(sequences)
        self.batch_worker.progress.connect(self.batch_progress.setValue)
        self.batch_worker.row_finished.connect(self.batch_row)
        self.batch_worker.finished.connect(lambda: self.batch_run.setEnabled(True))
        self.batch_worker.start()

    def batch_row(self, row, result, error):
        if error:
            self.table.setItem(row, 2, QTableWidgetItem("ERROR"))
            self.table.setItem(row, 3, QTableWidgetItem(error))
            return
        self.batch_results[row] = result
        values = [
            len(result.stems), len(result.qubo), result.prediction, result.reference,
            result.metrics["TP"], result.metrics["TN"], result.metrics["FP"],
            result.metrics["FN"], f"{result.metrics['MCC']:.4f}"
        ]
        for col, value in enumerate(values, 2):
            self.table.setItem(row, col, QTableWidgetItem(str(value)))

    def batch_detail(self, row, _column):
        if row >= len(self.batch_results) or self.batch_results[row] is None:
            return
        result = self.batch_results[row]
        self.tabs.setCurrentIndex(0)
        self.sequence.setPlainText(result.sequence)
        self.show_result(result)

    def export(self):
        valid = [r for r in self.batch_results if r is not None]
        if not valid:
            QMessageBox.warning(self, "No Results", "Run a batch prediction first.")
            return
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Results", "rna_batch_results.csv", "CSV Files (*.csv)"
        )
        if not filename:
            return
        with open(filename, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "Sequence", "Candidate Stems", "QUBO Coefficients", "QUBO Structure",
                "ViennaRNA Structure", "ViennaRNA MFE", "QUBO Energy",
                "TP", "TN", "FP", "FN", "MCC"
            ])
            for r in valid:
                w.writerow([
                    r.sequence, len(r.stems), len(r.qubo), r.prediction, r.reference,
                    r.mfe, r.solver_energy, r.metrics["TP"], r.metrics["TN"],
                    r.metrics["FP"], r.metrics["FN"], r.metrics["MCC"]
                ])
        QMessageBox.information(self, "Export Complete", f"Results saved to:\n{filename}")

    def apply_style(self):
        self.setStyleSheet("""
            QMainWindow, QWidget { background:#111318; color:#e8eaf0; font-family:"Noto Sans","Segoe UI",sans-serif; font-size:13px; }
            QGroupBox { border:1px solid #30343d; border-radius:10px; margin-top:12px; padding:14px; font-weight:bold; }
            QGroupBox::title { subcontrol-origin:margin; left:12px; padding:0 6px; color:#b9c0cc; }
            QTextEdit,QComboBox { background:#191c22; border:1px solid #363b45; border-radius:7px; padding:8px; color:#f0f2f6; }
            QPushButton { background:#2d6cdf; border:none; border-radius:7px; padding:9px 16px; color:white; font-weight:bold; }
            QPushButton:hover { background:#3d7bed; }
            QPushButton:disabled { background:#3a3e47; color:#858a94; }
            QTableWidget { background:#191c22; alternate-background-color:#15181d; gridline-color:#30343d; border:1px solid #30343d; }
            QHeaderView::section { background:#20242b; padding:8px; border:none; font-weight:bold; }
            #Header { font-size:27px; font-weight:bold; }
            #Subtitle { color:#9299a6; }
            #StatCard,#StructureCard { background:#191c22; border:1px solid #30343d; border-radius:9px; }
            #StatName { color:#9299a6; font-size:11px; }
            #StatValue { font-size:18px; font-weight:bold; }
            #CardTitle { font-size:15px; font-weight:bold; }
            #StructureText { background:#111318; border-radius:6px; padding:9px; font-family:monospace; font-size:15px; }
        """)


def run_gui():
    app = QApplication.instance() or QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(run_gui())
