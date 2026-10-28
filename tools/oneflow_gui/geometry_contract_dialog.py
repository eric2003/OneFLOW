"""Geometry contract explorer — table scenarios, not unit-test runner.

Uses tools/geometry_workbench (stub pipeline today; native C++ later).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QVBoxLayout,
)


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = (repository_root / "tools").resolve()
    tools_str = str(tools_root)
    if tools_str not in sys.path:
        sys.path.insert(0, tools_str)


class GeometryContractDialog(QDialog):
    """Browse and run table-level geometry contract scenarios inside Workbench."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        _ensure_tools_on_path(self._repository_root)

        self.setWindowTitle("Geometry Contracts")
        self.resize(900, 640)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Explore interface dual-map / split contracts as editable table scenarios. "
            "This is not the Unit Tests runner: it drives the same table pipeline the "
            "geometry refactor is carving out (stub now, native C++ later)."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        row = QHBoxLayout()
        row.addWidget(QLabel("Scenario"))
        self.scenario_combo = QComboBox()
        row.addWidget(self.scenario_combo, stretch=1)
        self.run_button = QPushButton("Run (stub)")
        self.run_button.clicked.connect(self._run_scenario)
        row.addWidget(self.run_button)
        root.addLayout(row)

        self.description = QLabel("")
        self.description.setWordWrap(True)
        root.addWidget(self.description)

        self.input_view = QPlainTextEdit()
        self.input_view.setReadOnly(True)
        self.input_view.setPlaceholderText("Scenario input tables (JSON)")
        root.addWidget(QLabel("Input"))
        root.addWidget(self.input_view, stretch=1)

        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        self.result_view.setPlaceholderText("Run result and topology sketch")
        root.addWidget(QLabel("Result"))
        root.addWidget(self.result_view, stretch=1)

        close_row = QHBoxLayout()
        close_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        close_row.addWidget(close_btn)
        root.addLayout(close_row)

        self._scenario_paths: list[Path] = []
        self._reload_scenarios()
        self.scenario_combo.currentIndexChanged.connect(self._on_scenario_changed)
        if self.scenario_combo.count():
            self._on_scenario_changed(0)

    def _reload_scenarios(self) -> None:
        from geometry_workbench.scenarios import list_scenario_paths

        self.scenario_combo.blockSignals(True)
        self.scenario_combo.clear()
        self._scenario_paths = list_scenario_paths(self._repository_root)
        for path in self._scenario_paths:
            self.scenario_combo.addItem(path.stem, str(path))
        self.scenario_combo.blockSignals(False)
        if not self._scenario_paths:
            self.description.setText(
                f"No scenarios under {self._repository_root / 'tools' / 'geometry_workbench' / 'scenarios'}"
            )

    def _current_scenario(self) -> dict[str, Any] | None:
        from geometry_workbench.scenarios import load_scenario

        index = self.scenario_combo.currentIndex()
        if index < 0 or index >= len(self._scenario_paths):
            return None
        return load_scenario(self._scenario_paths[index])

    def _on_scenario_changed(self, _index: int) -> None:
        scenario = self._current_scenario()
        if scenario is None:
            return
        self.description.setText(str(scenario.get("description") or scenario.get("id") or ""))
        self.input_view.setPlainText(json.dumps(scenario.get("input") or {}, indent=2))
        self.result_view.clear()

    def _run_scenario(self) -> None:
        from geometry_workbench.pipeline_stub import diff_tables, run_pipeline
        from geometry_workbench.visualize import render_ascii

        scenario = self._current_scenario()
        if scenario is None:
            QMessageBox.warning(self, "Geometry Contracts", "No scenario selected.")
            return
        try:
            actual = run_pipeline(scenario, backend="stub")
        except Exception as exc:  # noqa: BLE001 — surface to operator
            QMessageBox.critical(self, "Geometry Contracts", str(exc))
            return

        expected = scenario.get("expected") or {}
        mismatches = diff_tables(actual, expected)
        orphans = expected.get("orphanParentGlobalIds") or []
        lines = [
            f"scenario: {scenario.get('id', '')}",
            "backend: stub",
            "RESULT: PASS" if not mismatches else "RESULT: FAIL",
        ]
        for line in mismatches:
            lines.append(f"  - {line}")
        lines.append("")
        lines.append(render_ascii(actual, orphan_ids=orphans))
        notes = scenario.get("notes") or []
        if notes:
            lines.append("")
            lines.append("notes:")
            for note in notes:
                lines.append(f"  - {note}")
        self.result_view.setPlainText("\n".join(lines))
