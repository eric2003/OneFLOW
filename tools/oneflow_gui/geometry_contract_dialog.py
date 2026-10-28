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
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QInputDialog,
    QPushButton,
    QPlainTextEdit,
    QListWidget,
    QSplitter,
    QTabWidget,
    QWidget,
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
        self.resize(1120, 720)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Select a named contract scenario on the left. Edit input and expected "
            "tables in separate tabs; the result tab shows comparison output. "
            "Fixture mode does not execute production C++."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        command_row = QHBoxLayout()
        self.import_button = QPushButton("Create from GoogleTest...")
        self.import_button.clicked.connect(self._create_scenario_from_gtest)
        command_row.addWidget(self.import_button)
        self.save_button = QPushButton("Save scenario")
        self.save_button.clicked.connect(self._save_scenario)
        command_row.addWidget(self.save_button)
        self.run_button = QPushButton("Run scenario (fixture)")
        self.run_button.clicked.connect(self._run_scenario)
        command_row.addWidget(self.run_button)
        self.export_button = QPushButton("Export PNG...")
        self.export_button.setToolTip(
            "Optional static figure (matplotlib Agg); no interactive window"
        )
        self.export_button.clicked.connect(self._export_png)
        self.export_button.setEnabled(False)
        command_row.addWidget(self.export_button)
        command_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        command_row.addWidget(close_btn)
        root.addLayout(command_row)

        splitter = QSplitter()
        self.scenario_list = QListWidget()
        self.scenario_list.setMinimumWidth(220)
        self.scenario_list.currentRowChanged.connect(self._on_scenario_changed)
        splitter.addWidget(self.scenario_list)

        editor_panel = QWidget()
        editor_layout = QVBoxLayout(editor_panel)
        editor_layout.setContentsMargins(0, 0, 0, 0)
        self.description = QLabel("")
        self.description.setWordWrap(True)
        editor_layout.addWidget(self.description)
        self.tabs = QTabWidget()

        input_tab = QWidget()
        input_layout = QVBoxLayout(input_tab)
        input_header = QHBoxLayout()
        input_header.addWidget(QLabel("Input tables (editable JSON)"))
        input_header.addStretch(1)
        self.reset_input_button = QPushButton("Restore input")
        self.reset_input_button.clicked.connect(self._restore_input)
        input_header.addWidget(self.reset_input_button)
        input_layout.addLayout(input_header)
        self.input_view = QPlainTextEdit()
        self.input_view.setPlaceholderText("Edit scenario input tables (JSON)")
        input_layout.addWidget(self.input_view)
        self.tabs.addTab(input_tab, "Input")

        expected_tab = QWidget()
        expected_layout = QVBoxLayout(expected_tab)
        expected_header = QHBoxLayout()
        expected_header.addWidget(QLabel("Expected result (editable JSON)"))
        expected_header.addStretch(1)
        self.reset_expected_button = QPushButton("Restore expected")
        self.reset_expected_button.clicked.connect(self._restore_expected)
        expected_header.addWidget(self.reset_expected_button)
        expected_layout.addLayout(expected_header)
        self.expected_view = QPlainTextEdit()
        self.expected_view.setPlaceholderText("Edit expected result tables (JSON)")
        expected_layout.addWidget(self.expected_view)
        self.tabs.addTab(expected_tab, "Expected")

        result_tab = QWidget()
        result_layout = QVBoxLayout(result_tab)
        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        self.result_view.setPlaceholderText("Actual result, comparison, and notes")
        result_layout.addWidget(self.result_view)
        self.tabs.addTab(result_tab, "Result")

        editor_layout.addWidget(self.tabs)
        splitter.addWidget(editor_panel)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 4)
        root.addWidget(splitter, stretch=1)

        self._reload_scenarios()
        if self.scenario_list.count():
            self.scenario_list.setCurrentRow(0)

    def _reload_scenarios(self) -> None:
        from geometry_workbench.scenarios import list_scenario_paths

        self.scenario_list.blockSignals(True)
        self.scenario_list.clear()
        self._scenario_paths = list_scenario_paths(self._repository_root)
        from geometry_workbench.scenarios import load_scenario
        for path in self._scenario_paths:
            scenario = load_scenario(path)
            self.scenario_list.addItem(str(scenario.get("title") or path.stem))
        self.scenario_list.blockSignals(False)
        if not self._scenario_paths:
            self.description.setText(
                f"No scenarios under "
                f"{self._repository_root / 'tools' / 'geometry_workbench' / 'scenarios'}"
            )

    def _create_scenario_from_gtest(self) -> None:
        from geometry_workbench.scenario_import import write_scenario_skeleton
        from oneflow_gui.test_inventory import discover_tests

        records = [
            record for record in discover_tests(self._repository_root)
            if record.get("framework") == "GoogleTest" and not record.get("scenarioId")
        ]
        if not records:
            QMessageBox.information(
                self, "Geometry Contracts", "No unlinked GoogleTest declarations were discovered."
            )
            return

        labels = [
            f"{record['suite']}.{record['name']} — {record['source']}:{record['line']}"
            for record in records
        ]
        selected, accepted = QInputDialog.getItem(
            self,
            "Create scenario from GoogleTest",
            "Select a GoogleTest declaration. This creates a scenario skeleton; "
            "input and expected data must be completed or extracted separately.",
            labels,
            0,
            False,
        )
        if not accepted:
            return
        record = records[labels.index(selected)]
        scenario_dir = (
            self._repository_root / "tools" / "geometry_workbench" / "scenarios"
        )
        try:
            path = write_scenario_skeleton(record, scenario_dir)
        except (OSError, ValueError, FileExistsError) as exc:
            QMessageBox.warning(self, "Geometry Contracts", str(exc))
            return

        self._reload_scenarios()
        for index, scenario_path in enumerate(self._scenario_paths):
            if scenario_path.resolve() == path.resolve():
                self.scenario_list.setCurrentRow(index)
                break
        QMessageBox.information(
            self,
            "Geometry Contracts",
            "Created a traceable scenario skeleton. The input and expected tables "
            "are intentionally empty; this action does not execute C++ or infer "
            f"assertions.\n\n{path}",
        )

    def _current_scenario(self) -> dict[str, Any] | None:
        from geometry_workbench.scenarios import load_scenario

        index = self.scenario_list.currentRow()
        if index < 0 or index >= len(self._scenario_paths):
            return None
        return load_scenario(self._scenario_paths[index])

    def _on_scenario_changed(self, _index: int) -> None:
        scenario = self._current_scenario()
        if scenario is None:
            return
        description = str(scenario.get("description") or scenario.get("id") or "")
        linked_tests = [
            scenario.get("test") or {},
            scenario.get("integrationTest") or {},
        ]
        for test in linked_tests:
            if test:
                description += (
                    f"\nC++ regression (linked, not run here): "
                    f"{test.get('suite', '')}.{test.get('name', '')} "
                    f"({test.get('source', '')})"
                )
        self.description.setText(description)
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

    def _save_scenario(self) -> None:
        from geometry_workbench.scenarios import save_scenario_contract

        scenario = self._current_scenario()
        index = self.scenario_list.currentRow()
        if scenario is None or index < 0:
            QMessageBox.warning(self, "Geometry Contracts", "No scenario selected.")
            return
        try:
            input_tables = json.loads(self.input_view.toPlainText())
            expected_tables = json.loads(self.expected_view.toPlainText())
            save_scenario_contract(
                self._scenario_paths[index], scenario, input_tables, expected_tables
            )
        except (json.JSONDecodeError, OSError, ValueError) as exc:
            QMessageBox.warning(
                self, "Geometry Contracts", f"Could not save scenario contract: {exc}"
            )
            return
        QMessageBox.information(
            self,
            "Geometry Contracts",
            "Saved input and expected tables. Execution result was not changed.",
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
        has_topology_map = any(key in actual for key in ("l2g", "gI2Zid", "g2l"))
        self.export_button.setEnabled(has_topology_map)

        lines = [
            f"scenario: {scenario.get('id', '')}",
            "linked C++ regressions (not executed by this fixture run):",
        ]
        for test in (scenario.get("test") or {}, scenario.get("integrationTest") or {}):
            if test:
                lines.append(
                    f"  - {test.get('suite', '')}.{test.get('name', '')}"
                )
        lines.extend([
            "backend: fixture (does not execute the production C++ test)",
            "RESULT: PASS" if not mismatches else "RESULT: FAIL",
            (
                "Compared computed scenario output against the editable expected tables."
                if scenario.get("id") == "registered_child_resolution"
                else "Compared fixture output against the editable expected tables."
            ),
            (
                "The Python contract runner computed output from input; it does not execute production C++."
                if scenario.get("id") == "registered_child_resolution"
                else "Input JSON was parsed; this fixture backend does not recompute output from input."
            ),
        ])
        for line in mismatches:
            lines.append(f"  - {line}")
        lines.append("")
        if "l2g" in actual or "gI2Zid" in actual or "g2l" in actual:
            lines.append(render_ascii(actual, orphan_ids=orphans))
        else:
            lines.append("Computed result tables:")
            lines.append(json.dumps(actual, indent=2, ensure_ascii=True))
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
