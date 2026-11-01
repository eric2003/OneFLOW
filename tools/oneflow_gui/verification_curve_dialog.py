"""Theory-curve browser for intermediate numerical verification."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = (repository_root / "tools").resolve()
    tools_str = str(tools_root)
    if tools_str not in sys.path:
        sys.path.insert(0, tools_str)


class VerificationCurveDialog(QDialog):
    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        _ensure_tools_on_path(self._repository_root)

        from verification.curve_catalog import list_theory_curves

        self._curves = list_theory_curves()

        self.setWindowTitle("Verification Curves (Contract Workbench)")
        self.resize(900, 560)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Intermediate theory curves (closed-form). These are not full CFD "
            "cases. Shared product name: Contract Workbench "
            "(geometry table contracts + numerical theory curves)."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        body = QHBoxLayout()
        self.list = QListWidget()
        for c in self._curves:
            self.list.addItem(f"[{c.domain}] {c.title}")
        self.list.currentRowChanged.connect(self._on_select)
        body.addWidget(self.list, 1)

        right = QVBoxLayout()
        self.detail = QPlainTextEdit()
        self.detail.setReadOnly(True)
        right.addWidget(self.detail, 1)
        row = QHBoxLayout()
        self.export_button = QPushButton("Export PNG...")
        self.export_button.clicked.connect(self._export_png)
        row.addWidget(self.export_button)
        row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        row.addWidget(close_btn)
        right.addLayout(row)
        body.addLayout(right, 2)
        root.addLayout(body)

        if self._curves:
            self.list.setCurrentRow(0)

    def _on_select(self, row: int) -> None:
        if row < 0 or row >= len(self._curves):
            self.detail.clear()
            return
        c = self._curves[row]
        from verification.curve_catalog import csv_path

        path = csv_path(self._repository_root, c)
        self.detail.setPlainText(
            f"id: {c.id}\n"
            f"domain: {c.domain}\n"
            f"title: {c.title}\n\n"
            f"{c.description}\n\n"
            f"CSV: {path}\n"
            f"exists: {path.is_file()}\n\n"
            f"Proves: {c.proves}\n"
            f"Does not prove: {c.does_not_prove}\n"
        )

    def _export_png(self) -> None:
        row = self.list.currentRow()
        if row < 0:
            return
        c = self._curves[row]
        out, _ = QFileDialog.getSaveFileName(
            self, "Export theory curve PNG", f"{c.id}.png", "PNG (*.png)"
        )
        if not out:
            return
        try:
            from verification.curve_catalog import render_png

            render_png(self._repository_root, c, Path(out))
            QMessageBox.information(self, "Export", f"Wrote\n{out}")
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
