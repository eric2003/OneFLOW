"""Geometry contract explorer - table scenarios, not unit-test runner.

Uses tools/geometry_workbench (stub pipeline today; native C++ later).
Graphics export is optional and non-interactive (static PNG).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
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
        self._scenario_paths: list[Path] = []
        self._last_actual: dict[str, Any] | None = None
        self._last_orphans: list[int] = []

        self.setWindowTitle("Geometry Contracts")
        self.resize(900, 640)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Explore interface dual-map / split contracts as table scenarios. "
            "This is not the Unit Tests runner. Text PASS/FAIL is the default; "
            "PNG export is optional and does not open an interactive plot window."
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
        self.export_button = QPushButton("Export PNG...")
        self.export_button.setToolTip(
            "Optional static figure (matplotlib Agg); no interactive window"
        )
        self.export_button.clicked.connect(self._export_png)
        self.export_button.setEnabled(False)
        row.addWidget(self.export_button)
        root.addLayout(row)

        self.description = QLabel("")
        self.description.setWordWrap(True)
        root.addWidget(self.description)

        root.addWidget(QLabel("Input"))
        input_row = QHBoxLayout()
        input_row.addWidget(QLabel("Input tables (editable JSON)"))
        input_row.addStretch(1)
        self.reset_input_button = QPushButton("Restore scenario input")
        self.reset_input_button.clicked.connect(self._restore_input)
        input_row.addWidget(self.reset_input_button)
        root.addLayout(input_row)
        self.input_view = QPlainTextEdit()
        self.input_view.setPlaceholderText("Edit scenario input tables (JSON)")
        root.addWidget(self.input_view, stretch=1)

        expected_row = QHBoxLayout()
        expected_row.addWidget(QLabel("Expected tables (editable JSON)"))
        expected_row.addStretch(1)
        self.reset_expected_button = QPushButton("Restore scenario expected")
        self.reset_expected_button.clicked.connect(self._restore_expected)
        expected_row.addWidget(self.reset_expected_button)
        root.addLayout(expected_row)
        self.expected_view = QPlainTextEdit()
        self.expected_view.setPlaceholderText("Edit expected tables to probe the contract")
        root.addWidget(self.expected_view, stretch=1)

        root.addWidget(QLabel("Result"))
        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        self.result_view.setPlaceholderText("Mock pipeline result, comparison, and topology sketch")
        root.addWidget(self.result_view, stretch=1)

        close_row = QHBoxLayout()
        close_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        close_row.addWidget(close_btn)
        root.addLayout(close_row)

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
                f"No scenarios under "
                f"{self._repository_root / 'tools' / 'geometry_workbench' / 'scenarios'}"
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
        self.expected_view.setPlainText(
            json.dumps(scenario.get("expected") or {}, indent=2)
        )
        self.result_view.clear()
        self._last_actual = None
        self._last_orphans = []
        self.export_button.setEnabled(False)

    def _restore_input(self) -> None:
        scenario = self._current_scenario()
        if scenario is None:
            return
        self.input_view.setPlainText(json.dumps(scenario.get("input") or {}, indent=2))

    def _restore_expected(self) -> None:
        scenario = self._current_scenario()
        if scenario is None:
            return
        self.expected_view.setPlainText(
            json.dumps(scenario.get("expected") or {}, indent=2)
        )

    def _run_scenario(self) -> None:
        from geometry_workbench.pipeline_stub import diff_tables, run_pipeline
        from geometry_workbench.visualize import render_ascii

        scenario = self._current_scenario()
        if scenario is None:
            QMessageBox.warning(self, "Geometry Contracts", "No scenario selected.")
            return
        try:
            edited_input = json.loads(self.input_view.toPlainText())
            if not isinstance(edited_input, dict):
                raise ValueError("input tables must be a JSON object")
            scenario["input"] = edited_input
        except (json.JSONDecodeError, ValueError) as exc:
            QMessageBox.warning(
                self, "Geometry Contracts", f"Invalid input JSON: {exc}"
            )
            return

        try:
            actual = run_pipeline(scenario, backend="stub")
        except Exception as exc:  # noqa: BLE001 - surface to operator
            QMessageBox.critical(self, "Geometry Contracts", str(exc))
            return

        try:
            expected = json.loads(self.expected_view.toPlainText())
            if not isinstance(expected, dict):
                raise ValueError("expected tables must be a JSON object")
        except (json.JSONDecodeError, ValueError) as exc:
            QMessageBox.warning(
                self, "Geometry Contracts", f"Invalid expected JSON: {exc}"
            )
            return
        mismatches = diff_tables(actual, expected)
        orphans = actual.get("orphanParentGlobalIds") or []
        self._last_actual = actual
        self._last_orphans = list(orphans)
        self.export_button.setEnabled(True)

        lines = [
            f"scenario: {scenario.get('id', '')}",
            "backend: stub (fixture-backed, not the production geometry algorithm)",
            "RESULT: PASS" if not mismatches else "RESULT: FAIL",
            "Compared mock output against the editable expected tables.",
            "Input JSON was parsed; this fixture backend does not recompute output from input.",
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

    def _export_png(self) -> None:
        if not self._last_actual:
            QMessageBox.information(self, "Geometry Contracts", "Run a scenario first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export topology PNG", "geometry_contract.png", "PNG (*.png)"
        )
        if not path:
            return
        try:
            from geometry_workbench.visualize import save_topology_png

            out = save_topology_png(
                self._last_actual, Path(path), orphan_ids=self._last_orphans
            )
        except ImportError:
            QMessageBox.warning(
                self,
                "Geometry Contracts",
                "matplotlib is not installed. Text result is enough; PNG export is optional.",
            )
            return
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Geometry Contracts", str(exc))
            return
        QMessageBox.information(self, "Geometry Contracts", f"Wrote static figure:\n{out}")
