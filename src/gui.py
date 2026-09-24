import csv
import sys

import main
from vienna import Vienna

from PySide6.QtCore import Qt, QThread, Signal, QByteArray, QSize
from PySide6.QtGui import QImage, QPainter
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QApplication, QFileDialog, QFrame, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QMainWindow, QMessageBox, QProgressBar,
    QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QVBoxLayout, QWidget, QComboBox
)


class PredictionWorker(QThread):
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, sequence):
        super().__init__()
        self.sequence = sequence

    def run(self):
        try:
            # All scientific computation is delegated to main.py.
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

            self.progress.emit(
                int((i + 1) * 100 / len(self.sequences))
            )

        self.finished.emit()


class StructureCard(QFrame):
    """
    Displays a ViennaRNA SVG entirely from memory.

    No temporary image file is created. The SVG is kept as a string
    and passed directly to QSvgWidget.
    """

    def __init__(self, title, structure, svg_data):
        super().__init__()

        self.setObjectName("StructureCard")
        self.title = title
        self.structure = structure
        self.svg_data = svg_data

        layout = QVBoxLayout(self)

        title_label = QLabel(title)
        title_label.setObjectName("CardTitle")

        structure_label = QLabel(structure)
        structure_label.setObjectName("StructureText")
        structure_label.setAlignment(Qt.AlignCenter)
        structure_label.setTextInteractionFlags(
            Qt.TextSelectableByMouse
        )

        self.svg_widget = QSvgWidget()
        self.svg_widget.setMinimumSize(QSize(430, 300))
        self.svg_widget.setSizePolicy(
            self.svg_widget.sizePolicy().horizontalPolicy(),
            self.svg_widget.sizePolicy().verticalPolicy()
        )
        self.set_svg(svg_data)

        save_button = QPushButton("Save Image")
        save_button.clicked.connect(self.save_image)

        layout.addWidget(title_label)
        layout.addWidget(structure_label)
        layout.addWidget(self.svg_widget, 1)
        layout.addWidget(save_button)

    def set_svg(self, svg_data):
        if not isinstance(svg_data, str):
            raise TypeError(
                "Vienna.plot() must return SVG text, "
                f"got {type(svg_data).__name__}."
            )

        self.svg_data = svg_data.encode("utf-8")
        self.svg_widget.load(QByteArray(self.svg_data))

    def save_image(self):
        filename, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save RNA Structure Image",
            "",
            "SVG Image (*.svg);;PNG Image (*.png)"
        )

        if not filename:
            return

        try:
            if filename.lower().endswith(".png"):
                self.save_png(filename)
            else:
                if not filename.lower().endswith(".svg"):
                    filename += ".svg"

                with open(filename, "wb") as file:
                    file.write(self.svg_data)

            QMessageBox.information(
                self,
                "Image Saved",
                f"Structure image saved to:\n{filename}"
            )

        except Exception as exc:
            QMessageBox.critical(
                self,
                "Save Error",
                str(exc)
            )

    def save_png(self, filename):
        """
        Render the in-memory SVG directly into an in-memory QImage,
        then save that image to the user-selected path.
        """

        renderer = QSvgRenderer(QByteArray(self.svg_data))

        # Use a large enough raster canvas for a presentation-quality PNG.
        image = QImage(
            1600,
            1100,
            QImage.Format_ARGB32
        )
        image.fill(Qt.white)

        painter = QPainter(image)
        renderer.render(painter)
        painter.end()

        if not image.save(filename, "PNG"):
            raise RuntimeError(
                "Qt could not write the PNG file."
            )


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.batch_results = []

        self.setWindowTitle(
            "RNA Secondary Structure Predictor"
        )
        self.resize(1280, 850)

        self.tabs = QTabWidget()
        self.tabs.addTab(
            self.single_tab(),
            "Single Sequence"
        )
        self.tabs.addTab(
            self.batch_tab(),
            "Batch Analysis"
        )

        self.setCentralWidget(self.tabs)
        self.apply_style()

    # ---------------------------------------------------------
    # Single sequence
    # ---------------------------------------------------------

    def single_tab(self):

        page = QWidget()
        root = QVBoxLayout(page)

        header = QLabel(
            "RNA Secondary Structure Predictor"
        )
        header.setObjectName("Header")

        subtitle = QLabel(
            "QUBO-based folding compared with ViennaRNA MFE"
        )
        subtitle.setObjectName("Subtitle")

        root.addWidget(header)
        root.addWidget(subtitle)

        box = QGroupBox("RNA Sequence")
        layout = QVBoxLayout(box)

        self.sequence = QTextEdit()
        self.sequence.setPlaceholderText(
            "Enter RNA sequence, e.g. "
            "GGCGAAAUCGCCUUUGGCGAAAUCGCC"
        )
        self.sequence.setMaximumHeight(100)

        self.predict = QPushButton(
            "Predict Structure"
        )
        self.predict.clicked.connect(
            self.run_single
        )

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
        self.results_layout = QVBoxLayout(
            self.results
        )

        label = QLabel(
            "Enter an RNA sequence and click "
            "Predict Structure."
        )
        label.setAlignment(Qt.AlignCenter)

        self.results_layout.addWidget(label)

        self.scroll.setWidget(self.results)
        root.addWidget(self.scroll)

        return page

    def run_single(self):

        sequence = "".join(
            self.sequence.toPlainText().split()
        ).upper()

        if not sequence:
            QMessageBox.warning(
                self,
                "Missing sequence",
                "Please enter an RNA sequence."
            )
            return

        self.predict.setEnabled(False)
        self.progress.show()

        self.worker = PredictionWorker(sequence)

        self.worker.finished.connect(
            self.show_result
        )
        self.worker.error.connect(
            self.single_error
        )

        self.worker.start()

    def single_error(self, message):

        self.predict.setEnabled(True)
        self.progress.hide()

        QMessageBox.critical(
            self,
            "Prediction Error",
            message
        )

    def clear_results(self):

        while self.results_layout.count():

            item = self.results_layout.takeAt(0)

            if item.widget():
                item.widget().deleteLater()

    # ---------------------------------------------------------
    # Results
    # ---------------------------------------------------------

    def show_result(self, result):

        self.predict.setEnabled(True)
        self.progress.hide()

        self.clear_results()

        grid = QGridLayout()

        stats = [
            (
                "Sequence Length",
                len(result.sequence)
            ),
            (
                "Candidate Stems",
                len(result.stems)
            ),
            (
                "QUBO Coefficients",
                len(result.qubo)
            ),
            (
                "ViennaRNA MFE",
                f"{result.mfe:.2f} kcal/mol"
            ),
            (
                "QUBO Energy",
                f"{result.solver_energy:.4f}"
            ),
            (
                "TP",
                result.metrics["TP"]
            ),
            (
                "TN",
                result.metrics["TN"]
            ),
            (
                "FP",
                result.metrics["FP"]
            ),
            (
                "FN",
                result.metrics["FN"]
            ),
            (
                "MCC",
                f"{result.metrics['MCC']:.4f}"
            ),
        ]

        for i, (name, value) in enumerate(stats):

            card = QFrame()
            card.setObjectName("StatCard")

            layout = QVBoxLayout(card)

            name_label = QLabel(name)
            name_label.setObjectName(
                "StatName"
            )

            value_label = QLabel(
                str(value)
            )
            value_label.setObjectName(
                "StatValue"
            )

            layout.addWidget(name_label)
            layout.addWidget(value_label)

            grid.addWidget(
                card,
                i // 5,
                i % 5
            )

        self.results_layout.addLayout(grid)

        # -----------------------------------------------------
        # ViennaRNA-native structure drawings
        # -----------------------------------------------------

        vienna = Vienna(result.sequence)

        qubo_svg = vienna.plot(
            result.prediction
        )

        reference_svg = vienna.plot(
            result.reference
        )

        comparison = QGroupBox(
            "Structure Comparison"
        )

        row = QHBoxLayout(comparison)

        row.addWidget(
            StructureCard(
                "QUBO Predicted Structure",
                result.prediction,
                qubo_svg
            )
        )

        row.addWidget(
            StructureCard(
                "ViennaRNA MFE Structure",
                result.reference,
                reference_svg
            )
        )

        self.results_layout.addWidget(
            comparison
        )

        note = QLabel(
            "TP/TN/FP/FN are adjacency-matrix "
            "entry counts. MCC summarizes the "
            "prediction against the ViennaRNA reference."
        )

        note.setWordWrap(True)

        self.results_layout.addWidget(
            note
        )

        self.results_layout.addStretch()

    # ---------------------------------------------------------
    # Batch
    # ---------------------------------------------------------

    def batch_tab(self):

        page = QWidget()
        root = QVBoxLayout(page)

        header = QLabel(
            "Batch RNA Analysis"
        )
        header.setObjectName("Header")

        root.addWidget(header)

        root.addWidget(
            QLabel(
                "Upload a CSV or Excel file containing "
                "one RNA sequence per row."
            )
        )

        row = QHBoxLayout()

        browse = QPushButton(
            "Upload CSV / Excel"
        )
        browse.clicked.connect(
            self.load_file
        )

        self.file_label = QLabel(
            "No file selected"
        )

        self.columns = QComboBox()
        self.columns.setEnabled(False)

        self.batch_run = QPushButton(
            "Run Batch Prediction"
        )
        self.batch_run.setEnabled(False)
        self.batch_run.clicked.connect(
            self.run_batch
        )

        row.addWidget(browse)
        row.addWidget(
            self.file_label,
            1
        )
        row.addWidget(
            QLabel("Sequence column:")
        )
        row.addWidget(self.columns)
        row.addWidget(self.batch_run)

        root.addLayout(row)

        self.batch_progress = QProgressBar()
        root.addWidget(
            self.batch_progress
        )

        self.table = QTableWidget()

        self.table.setColumnCount(11)

        self.table.setHorizontalHeaderLabels([
            "ID",
            "Sequence",
            "Stems",
            "QUBO Coefficients",
            "QUBO Structure",
            "ViennaRNA",
            "TP",
            "TN",
            "FP",
            "FN",
            "MCC"
        ])

        self.table.setAlternatingRowColors(
            True
        )

        self.table.cellDoubleClicked.connect(
            self.batch_detail
        )

        root.addWidget(self.table)

        export = QPushButton(
            "Export Results as CSV"
        )
        export.clicked.connect(
            self.export
        )

        root.addWidget(export)

        self.batch_data = []

        return page

    def load_file(self):

        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Select RNA sequence file",
            "",
            "Data Files (*.csv *.xlsx *.xls)"
        )

        if not filename:
            return

        try:

            suffix = (
                filename.lower()
                .split(".")[-1]
            )

            if suffix == "csv":

                with open(
                    filename,
                    "r",
                    encoding="utf-8-sig",
                    newline=""
                ) as file:

                    reader = csv.DictReader(
                        file
                    )

                    headers = (
                        reader.fieldnames
                        or []
                    )

                    rows = list(reader)

            else:

                import pandas as pd

                dataframe = (
                    pd.read_excel(
                        filename
                    ).fillna("")
                )

                headers = [
                    str(column)
                    for column
                    in dataframe.columns
                ]

                rows = dataframe.to_dict(
                    "records"
                )

            self.batch_data = rows

            self.file_label.setText(
                filename
            )

            self.columns.clear()

            self.columns.addItems(
                headers
            )

            self.columns.setEnabled(
                True
            )

            self.batch_run.setEnabled(
                True
            )

        except Exception as exc:

            QMessageBox.critical(
                self,
                "File Error",
                str(exc)
            )

    def run_batch(self):

        column = (
            self.columns.currentText()
        )

        sequences = [
            "".join(
                str(
                    row.get(
                        column,
                        ""
                    )
                ).split()
            ).upper()
            for row
            in self.batch_data
        ]

        sequences = [
            sequence
            for sequence
            in sequences
            if sequence
        ]

        if not sequences:

            QMessageBox.warning(
                self,
                "No sequences",
                "No non-empty RNA sequences were found."
            )

            return

        self.batch_results = (
            [None] * len(sequences)
        )

        self.table.setRowCount(
            len(sequences)
        )

        for i, sequence in enumerate(
            sequences
        ):

            self.table.setItem(
                i,
                0,
                QTableWidgetItem(
                    str(i + 1)
                )
            )

            self.table.setItem(
                i,
                1,
                QTableWidgetItem(
                    sequence
                )
            )

        self.batch_run.setEnabled(
            False
        )

        self.batch_progress.setValue(
            0
        )

        self.batch_worker = BatchWorker(
            sequences
        )

        self.batch_worker.progress.connect(
            self.batch_progress.setValue
        )

        self.batch_worker.row_finished.connect(
            self.batch_row
        )

        self.batch_worker.finished.connect(
            lambda:
            self.batch_run.setEnabled(
                True
            )
        )

        self.batch_worker.start()

    def batch_row(
        self,
        row,
        result,
        error
    ):

        if error:

            self.table.setItem(
                row,
                2,
                QTableWidgetItem(
                    "ERROR"
                )
            )

            self.table.setItem(
                row,
                3,
                QTableWidgetItem(
                    error
                )
            )

            return

        self.batch_results[row] = result

        values = [
            len(result.stems),
            len(result.qubo),
            result.prediction,
            result.reference,
            result.metrics["TP"],
            result.metrics["TN"],
            result.metrics["FP"],
            result.metrics["FN"],
            f"{result.metrics['MCC']:.4f}"
        ]

        for column, value in enumerate(
            values,
            2
        ):

            self.table.setItem(
                row,
                column,
                QTableWidgetItem(
                    str(value)
                )
            )

    def batch_detail(
        self,
        row,
        _column
    ):

        if (
            row >= len(
                self.batch_results
            )
            or
            self.batch_results[row]
            is None
        ):
            return

        result = (
            self.batch_results[row]
        )

        self.tabs.setCurrentIndex(
            0
        )

        self.sequence.setPlainText(
            result.sequence
        )

        self.show_result(
            result
        )

    # ---------------------------------------------------------
    # CSV export
    # ---------------------------------------------------------

    def export(self):

        valid = [
            result
            for result
            in self.batch_results
            if result is not None
        ]

        if not valid:

            QMessageBox.warning(
                self,
                "No Results",
                "Run a batch prediction first."
            )

            return

        filename, _ = (
            QFileDialog.getSaveFileName(
                self,
                "Save Results",
                "rna_batch_results.csv",
                "CSV Files (*.csv)"
            )
        )

        if not filename:
            return

        with open(
            filename,
            "w",
            encoding="utf-8",
            newline=""
        ) as file:

            writer = csv.writer(
                file
            )

            writer.writerow([
                "Sequence",
                "Candidate Stems",
                "QUBO Coefficients",
                "QUBO Structure",
                "ViennaRNA Structure",
                "ViennaRNA MFE",
                "QUBO Energy",
                "TP",
                "TN",
                "FP",
                "FN",
                "MCC"
            ])

            for result in valid:

                writer.writerow([
                    result.sequence,
                    len(result.stems),
                    len(result.qubo),
                    result.prediction,
                    result.reference,
                    result.mfe,
                    result.solver_energy,
                    result.metrics["TP"],
                    result.metrics["TN"],
                    result.metrics["FP"],
                    result.metrics["FN"],
                    result.metrics["MCC"]
                ])

        QMessageBox.information(
            self,
            "Export Complete",
            f"Results saved to:\n{filename}"
        )

    # ---------------------------------------------------------
    # Styling
    # ---------------------------------------------------------

    def apply_style(self):

        self.setStyleSheet("""
            QMainWindow, QWidget {
                background: #111318;
                color: #e8eaf0;
                font-family:
                    "Noto Sans",
                    "Segoe UI",
                    sans-serif;
                font-size: 13px;
            }

            QGroupBox {
                border: 1px solid #30343d;
                border-radius: 10px;
                margin-top: 12px;
                padding: 14px;
                font-weight: bold;
            }

            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px;
                color: #b9c0cc;
            }

            QTextEdit,
            QComboBox {
                background: #191c22;
                border: 1px solid #363b45;
                border-radius: 7px;
                padding: 8px;
                color: #f0f2f6;
            }

            QPushButton {
                background: #2d6cdf;
                border: none;
                border-radius: 7px;
                padding: 9px 16px;
                color: white;
                font-weight: bold;
            }

            QPushButton:hover {
                background: #3d7bed;
            }

            QPushButton:disabled {
                background: #3a3e47;
                color: #858a94;
            }

            QTableWidget {
                background: #191c22;
                alternate-background-color: #15181d;
                gridline-color: #30343d;
                border: 1px solid #30343d;
            }

            QHeaderView::section {
                background: #20242b;
                padding: 8px;
                border: none;
                font-weight: bold;
            }

            #Header {
                font-size: 27px;
                font-weight: bold;
            }

            #Subtitle {
                color: #9299a6;
            }

            #StatCard,
            #StructureCard {
                background: #191c22;
                border: 1px solid #30343d;
                border-radius: 9px;
            }

            #StatName {
                color: #9299a6;
                font-size: 11px;
            }

            #StatValue {
                font-size: 18px;
                font-weight: bold;
            }

            #CardTitle {
                font-size: 15px;
                font-weight: bold;
            }

            #StructureText {
                background: #111318;
                border-radius: 6px;
                padding: 9px;
                font-family: monospace;
                font-size: 15px;
            }
        """)


def run_gui():

    app = (
        QApplication.instance()
        or QApplication(sys.argv)
    )

    window = MainWindow()
    window.show()

    return app.exec()


if __name__ == "__main__":
    raise SystemExit(
        run_gui()
    )
