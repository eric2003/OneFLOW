"""Browse source-discovered tests separately from scenario execution results."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QVBoxLayout,
)


class RepositoryTestInventoryDialog(QDialog):
    """Show the repository test catalog with explicit framework filtering."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        from .test_inventory import discover_tests

        self.setWindowTitle("Repository Test Inventory")
        self.resize(1050, 680)
        self._records = discover_tests(repository_root.resolve())
        self._visible_records: list[dict[str, Any]] = []

        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "This is a source catalog, not a run result. Input and expected "
                "values appear only for explicitly linked contracts; result is "
                "populated only by a real execution backend."
            )
        )

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Test framework"))
        self.framework_filter = QComboBox()
        self.framework_filter.addItem("All frameworks", "all")
        for framework in sorted({record["framework"] for record in self._records}):
            self.framework_filter.addItem(framework, framework)
        self.framework_filter.currentIndexChanged.connect(self._refresh_records)
        filter_row.addWidget(self.framework_filter)
        filter_row.addStretch(1)
        layout.addLayout(filter_row)

        splitter = QSplitter()
        self.test_list = QListWidget()
        self.test_list.setMinimumWidth(420)
        self.test_list.currentRowChanged.connect(self._show_record)
        splitter.addWidget(self.test_list)

        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setPlaceholderText(
            "Select a test to inspect input, expected, result, and source."
        )
        splitter.addWidget(self.details)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        layout.addWidget(splitter, 1)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        button_row.addWidget(close_button)
        layout.addLayout(button_row)

        self._refresh_records()

    def _refresh_records(self, _index: int = 0) -> None:
        framework = self.framework_filter.currentData()
        self._visible_records = [
            record
            for record in self._records
            if framework == "all" or record["framework"] == framework
        ]

        self.test_list.blockSignals(True)
        self.test_list.clear()
        self.test_list.addItem(
            f"{len(self._visible_records)} discovered tests; source scan only"
        )
        for record in self._visible_records:
            self.test_list.addItem(
                f"[{record['framework']}] {record['suite']}.{record['name']} "
                f"- {record['source']}:{record['line']}"
            )
        self.test_list.blockSignals(False)
        self.test_list.setCurrentRow(1 if self._visible_records else 0)

    def _show_record(self, row: int) -> None:
        index = row - 1
        if index < 0 or index >= len(self._visible_records):
            self.details.setPlainText(
                "Choose a test. GoogleTest and unittest are filterable separately. "
                "Discovery does not execute tests or invent missing contract data."
            )
            return

        self.details.setPlainText(
            json.dumps(
                {"schemaVersion": 1, **self._visible_records[index]},
                indent=2,
            )
        )
