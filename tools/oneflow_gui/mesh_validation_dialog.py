"""Dedicated Qt workbench for end-to-end Mesh-to-CFD validation."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtWidgets import (
    QComboBox, QDialog, QFormLayout, QHBoxLayout, QHeaderView, QLabel,
    QFileDialog, QMessageBox, QPlainTextEdit, QPushButton, QSpinBox, QSplitter,
    QTabWidget, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget,
)

from grid.stl_to_gmsh.run_report import write_run_report

from .executable_discovery import executable_candidates
from .mesh_validation import (
    cfd_template_choices, existing_mesh_case_paths, mesh_case_group_categories,
    mesh_case_lifecycle_summary, mesh_case_requires_oneflow_executable,
    run_report_stage_rows,
    geometry_lineage_rows,
    _collect_run_artifacts, _artifact_file_state,
    recover_mesh_validation_run, list_mesh_validation_runs,
)
from .process_control import process_group_options, terminate_process_tree
from .user_data_paths import UserDataPaths


class MeshValidationWorker(QThread):
    """Run all selected validation Cases, continuing after individual failures."""

    output_line = Signal(str)
    completed = Signal(bool, str)
    report_ready = Signal(str, object)

    def __init__(
        self, commands: list[tuple[str, tuple[str, ...], Path]],
        log_directory: Path,
    ) -> None:
        super().__init__()
        self._commands = commands
        self._log_directory = Path(log_directory)
        self._cancel = threading.Event()
        self._process: subprocess.Popen[str] | None = None
        self._lock = threading.Lock()

    def request_stop(self) -> None:
        self._cancel.set()
        with self._lock:
            process = self._process
        if process is not None:
            terminate_process_tree(process)

    def run(self) -> None:
        passed: list[str] = []
        failed: list[str] = []
        self._log_directory.mkdir(parents=True, exist_ok=True)
        for label, arguments, working_directory in self._commands:
            if self._cancel.is_set():
                break
            safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("._") or "case"
            log_path = self._log_directory / f"{safe_label}-{uuid.uuid4().hex[:12]}.log"
            case_output: list[str] = []
            with log_path.open("w", encoding="utf-8", newline="\n") as log_stream:
                def record(line: str) -> None:
                    log_stream.write(line + "\n")
                    log_stream.flush()
                    self.output_line.emit(line)

                record(f"=== {label} ===")
                record("Started (UTC): " + datetime.now(timezone.utc).isoformat())
                record("Working directory: " + str(working_directory))
                command_text = subprocess.list2cmdline(list(arguments))
                record("Command: " + command_text)
                try:
                    with subprocess.Popen(
                        list(arguments), cwd=working_directory,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                        text=True, encoding="utf-8", errors="replace", bufsize=1,
                        **process_group_options(),
                    ) as process:
                        with self._lock:
                            self._process = process
                        try:
                            if process.stdout is not None:
                                for line in process.stdout:
                                    rendered_line = line.rstrip("\r\n")
                                    case_output.append(rendered_line)
                                    record(rendered_line)
                                    if self._cancel.is_set():
                                        terminate_process_tree(process)
                                        break
                            return_code = process.wait()
                        finally:
                            with self._lock:
                                if self._process is process:
                                    self._process = None
                except OSError as error:
                    record(f"ERROR: {error}")
                    self.output_line.emit(f"Full execution log: {log_path}")
                    failed.append(label)
                    continue
                if self._cancel.is_set():
                    record("Run cancelled before report finalization.")
                    break
                report_path = next(
                    (
                        Path(line.strip().split(": ", 1)[1]).expanduser()
                        for line in reversed(case_output)
                        if line.strip().startswith("report: ")
                    ),
                    None,
                )
                report = None
                if report_path is not None and report_path.is_file():
                    try:
                        report = json.loads(report_path.read_text(encoding="utf-8"))
                    except (OSError, json.JSONDecodeError) as error:
                        record(f"WARNING: cannot read run report: {error}")
                outcome = "PASSED" if return_code == 0 else f"FAILED (exit {return_code})"
                record(f"{outcome}: {label}")
                record("Finished (UTC): " + datetime.now(timezone.utc).isoformat())
            if isinstance(report, dict) and report.get("schema") == "oneflow.run-report":
                durable_log = log_path
                if report_path is not None:
                    durable_log = report_path.parent / "mesh-validation-console.log"
                    try:
                        durable_log.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(log_path, durable_log)
                        log_path.unlink(missing_ok=True)
                    except OSError:
                        durable_log = log_path
                console_reference = str(durable_log)
                if report_path is not None:
                    try:
                        console_reference = durable_log.relative_to(report_path.parent).as_posix()
                    except ValueError:
                        pass
                execution = report.setdefault("execution", {})
                execution.update({
                    "command": list(arguments),
                    "working_directory": str(working_directory),
                    "console_log": console_reference,
                    "return_code": return_code,
                })
                report.setdefault("artifacts", {})["console_log"] = console_reference
                if report_path is not None:
                    report["artifacts"]["files"] = _collect_run_artifacts(
                        report_path.parent, report_path
                    )
                if report_path is not None:
                    try:
                        write_run_report(report_path, report)
                    except OSError as error:
                        self.output_line.emit(f"WARNING: cannot update report with console log: {error}")
                self.report_ready.emit(label, report)
            else:
                self.output_line.emit(
                    f"WARNING: no structured run report was found; full log: {log_path}"
                )
            if return_code == 0:
                passed.append(label)
            else:
                failed.append(label)
        cancelled = self._cancel.is_set()
        summary = [
            "Mesh-to-CFD Validation summary:",
            f"  Selected: {len(self._commands)}",
            f"  Passed: {len(passed)}",
            f"  Failed: {len(failed)}",
            f"  Not run: {len(self._commands) - len(passed) - len(failed)}",
        ]
        if failed:
            summary.append("  Failed Cases: " + ", ".join(failed))
        if cancelled:
            summary.append("  Run cancelled.")
        self.completed.emit(not failed and not cancelled, "\n".join(summary))

class MeshValidationDialog(QDialog):
    """Select mesh Cases and run their declared validation milestone."""

    def __init__(
        self, repository_root: Path, *, executable: str = "",
        mesh_file: Path | None = None, project_directory: Path | None = None,
        mesh_case_id: str = "stl_to_gmsh_sphere_automatic_tetrahedra",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self._worker: MeshValidationWorker | None = None
        self._mesh_file = Path(mesh_file).expanduser().resolve() if mesh_file else None
        self._project_directory = (
            Path(project_directory).expanduser().resolve()
            if project_directory else (self._mesh_file.parent if self._mesh_file else None)
        )
        self._mesh_case_id = mesh_case_id
        self.setWindowTitle("Mesh-to-CFD Validation")
        self.resize(1180, 760)
        self.setMinimumSize(900, 620)

        layout = QVBoxLayout(self)
        introduction = QLabel(
            "Validate the current MSH file without regenerating it. The mesh is "
            "copied into an isolated run directory, then passed through CGNS "
            "conversion and OneFLOW ingestion. The CFD stage is a short runtime "
            "smoke test, not a numerical-accuracy claim."
            if self._mesh_file is not None else
            "Select a native C++ or Python STL-to-Gmsh Case, or an existing "
            "Gmsh mesh produced in the Mesh Generation workbench. Cases with an "
            "incremental lifecycle contract show their current verification "
            "milestone beside the name; hover for implemented and planned stages. "
            "Only stages in the declared milestone are executed. Cases without a "
            "lifecycle contract retain the complete legacy workflow. Existing "
            "meshes are listed separately under Existing Gmsh Mesh; for 3D meshes, "
            "the consumer dimension is adjusted. A CFD smoke test is not a claim "
            "of numerical accuracy. Re-running a Case clears only its previous "
            "validation directory."
        )
        introduction.setWordWrap(True)
        layout.addWidget(introduction)

        self.current_mesh_label = QLabel(
            "Current mesh: " + str(self._mesh_file) if self._mesh_file is not None else ""
        )
        self.current_mesh_label.setWordWrap(True)
        if self._mesh_file is not None:
            layout.addWidget(self.current_mesh_label)

        self.case_tree = QTreeWidget()
        # Selection, tree expansion and case names belong to one column. A
        # separate empty "Run" column wastes space and misleads users into
        # thinking that it is a run-status indicator.
        self.case_tree.setColumnCount(1)
        self.case_tree.setHeaderLabels(["STL / Mesh Case"])
        self.case_tree.setTreePosition(0)
        self.case_tree.setRootIsDecorated(True)
        self.case_tree.setAlternatingRowColors(True)
        case_header = self.case_tree.header()
        case_header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        case_header.setStretchLastSection(True)
        case_header.setMinimumSectionSize(240)
        self._case_items: dict[str, QTreeWidgetItem] = {}
        try:
            dimension_groups = mesh_case_group_categories(self._repository_root)
            for dimension_label, categories in dimension_groups:
                dimension_item = QTreeWidgetItem(self.case_tree)
                dimension_item.setText(0, dimension_label)
                dimension_item.setToolTip(0, dimension_label)
                dimension_item.setFlags(
                    dimension_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                )
                for category_label, choices in categories:
                    category_item = QTreeWidgetItem(dimension_item)
                    category_item.setText(0, category_label)
                    category_item.setToolTip(0, category_label)
                    category_item.setFlags(
                        category_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                    )
                    for case_id, title in choices:
                        item = QTreeWidgetItem(category_item)
                        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                        lifecycle = mesh_case_lifecycle_summary(
                            self._repository_root, case_id
                        )
                        if lifecycle is None:
                            item.setText(0, title)
                            item.setToolTip(0, title)
                        else:
                            scope = lifecycle["verification_scope"]
                            item.setText(0, f"{title}  [Milestone: {scope}]")
                            stage_lines = []
                            for stage in lifecycle["stages"]:
                                scope_label = (
                                    "in current scope"
                                    if stage["in_current_scope"] else "future stage"
                                )
                                variants = ", ".join(
                                    f'{variant["id"]}={variant["status"]}'
                                    for variant in stage["implementations"]
                                ) or "no implementation variant"
                                stage_lines.append(
                                    f'{stage["id"]} ({stage["type"]}): '
                                    f'{stage["implementation_status"]} / '
                                    f'{stage["verification_status"]}; '
                                    f'{scope_label}; {variants}'
                                )
                            item.setToolTip(
                                0,
                                title + "\nCurrent verification milestone: " + scope
                                + "\n\nLifecycle stages:\n"
                                + "\n".join(stage_lines),
                            )
                        item.setData(0, Qt.ItemDataRole.UserRole, case_id)
                        item.setCheckState(0, Qt.CheckState.Checked)
                        self._case_items[case_id] = item
                    category_item.setCheckState(0, Qt.CheckState.Checked)
                dimension_item.setCheckState(0, Qt.CheckState.Checked)
            self.case_tree.itemChanged.connect(self._update_group_state)
            self.case_tree.expandToDepth(1)
        except (OSError, ValueError, KeyError) as error:
            self.case_tree.setEnabled(False)
            layout.addWidget(QLabel(f"Mesh Case catalog unavailable: {error}"))
        self.case_tree.setMinimumHeight(180)
        details_panel = QWidget(self)
        content_layout = QVBoxLayout(details_panel)
        content_layout.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Vertical, self)
        splitter.setChildrenCollapsible(False)
        splitter.setHandleWidth(7)
        splitter.addWidget(self.case_tree)
        splitter.addWidget(details_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([330, 370])
        layout.addWidget(splitter, 1)
        self.case_tree.setVisible(self._mesh_file is None)

        selection_actions = QHBoxLayout()
        select_all = QPushButton("Select all")
        select_all.clicked.connect(lambda: self._set_all(True))
        select_none = QPushButton("Select none")
        select_none.clicked.connect(lambda: self._set_all(False))
        for selection_button in (select_all, select_none):
            selection_button.setFixedHeight(26)
            selection_button.setStyleSheet("padding: 2px 8px;")
        selection_actions.addWidget(select_all)
        selection_actions.addWidget(select_none)
        selection_actions.addStretch(1)
        content_layout.addLayout(selection_actions)
        select_all.setVisible(self._mesh_file is None)
        select_none.setVisible(self._mesh_file is None)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.executable_combo.setMinimumContentsLength(28)
        self.executable_combo.setToolTip(
            "Path to the OneFLOW executable. The field can be resized with the dialog."
        )
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if executable:
            self.executable_combo.setEditText(executable)
        elif self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        self.template_combo = QComboBox()
        for template_id, title in cfd_template_choices():
            self.template_combo.addItem(title, template_id)
        self.steps_spin = QSpinBox()
        self.steps_spin.setRange(1, 100000)
        self.steps_spin.setValue(2)
        self.steps_spin.setToolTip("Short CFD solver smoke-test iteration count.")
        self.executable_combo.setPlaceholderText(
            "Required only when selected milestones include OneFLOW stages"
        )
        self.executable_combo.setToolTip(
            "A Python geometry/Gmsh-only milestone does not need OneFLOW. "
            "Select an executable when a selected Case includes conversion, "
            "ingestion, or CFD execution."
        )
        form.addRow("OneFLOW executable (if needed)", self.executable_combo)
        template_label = QLabel("CFD consumer template")
        template_label.setToolTip(
            "For 3D meshes, the runtime dimension is overridden from mesh metadata."
        )
        form.addRow(template_label, self.template_combo)
        form.addRow("CFD smoke-test steps", self.steps_spin)
        content_layout.addLayout(form)

        workspace = (
            self._project_directory / "validation"
            if self._mesh_file is not None and self._project_directory is not None
            else UserDataPaths.discover().mesh_validation
        )
        self.workspace_label = QLabel("Validation output: " + str(workspace))
        self.workspace_label.setWordWrap(True)
        content_layout.addWidget(self.workspace_label)

        self.result_tabs = QTabWidget()
        summary_page = QWidget(self)
        summary_layout = QVBoxLayout(summary_page)
        summary_layout.setContentsMargins(4, 4, 4, 4)
        self.result_summary = QLabel("No validation results yet.")
        self.result_summary.setWordWrap(True)
        summary_layout.addWidget(self.result_summary)
        self.result_tree = QTreeWidget()
        self.result_tree.setColumnCount(3)
        self.result_tree.setHeaderLabels(["Case / stage", "Status", "Details"])
        self.result_tree.setRootIsDecorated(True)
        self.result_tree.setAlternatingRowColors(True)
        result_header = self.result_tree.header()
        result_header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        result_header.setSectionResizeMode(1, QHeaderView.ResizeMode.Interactive)
        result_header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        result_header.setStretchLastSection(False)
        result_header.resizeSection(0, 220)
        result_header.resizeSection(1, 100)
        summary_layout.addWidget(self.result_tree, 1)
        self.result_tabs.addTab(summary_page, "Results")

        self.artifact_tree = QTreeWidget()
        self.artifact_tree.setColumnCount(3)
        self.artifact_tree.setHeaderLabels(["Artifact", "Availability", "Path / Details"])
        self.artifact_tree.setAlternatingRowColors(True)
        self.artifact_tree.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        self.artifact_tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.artifact_tree.header().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.result_tabs.addTab(self.artifact_tree, "Artifacts")
        self._reported_outcomes: list[str] = []

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.log_view.setPlaceholderText("Validation output will appear here.")
        self.result_tabs.addTab(self.log_view, "Console")

        self.raw_report_view = QPlainTextEdit()
        self.raw_report_view.setReadOnly(True)
        self.raw_report_view.setPlaceholderText("Open a validation report to inspect all recorded fields.")
        self.result_tabs.addTab(self.raw_report_view, "Raw Report")

        history_page = QWidget(self)
        history_layout = QVBoxLayout(history_page)
        history_actions = QHBoxLayout()
        history_actions.addWidget(QLabel("Previous validation runs"))
        history_actions.addStretch(1)
        self.refresh_history_button = QPushButton("Refresh")
        self.refresh_history_button.clicked.connect(self._refresh_run_history)
        history_actions.addWidget(self.refresh_history_button)
        history_layout.addLayout(history_actions)
        self.history_tree = QTreeWidget()
        self.history_tree.setColumnCount(4)
        self.history_tree.setHeaderLabels(["Mesh Case", "Status", "Started (UTC)", "Run directory"])
        self.history_tree.setAlternatingRowColors(True)
        self.history_tree.setRootIsDecorated(False)
        self.history_tree.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        self.history_tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.history_tree.header().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.history_tree.header().setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.history_tree.itemDoubleClicked.connect(self._open_selected_history_run)
        history_layout.addWidget(self.history_tree, 1)
        self.result_tabs.addTab(history_page, "Run History")
        content_layout.addWidget(self.result_tabs, 1)
        self._history_workspace = workspace
        self._refresh_run_history()

        buttons = QHBoxLayout()
        self.open_run_button = QPushButton("Open Run...")
        self.open_run_button.setToolTip("Reconstruct the result views from an existing mesh-validation-report.json.")
        self.open_run_button.clicked.connect(self._open_run_report)
        buttons.addWidget(self.open_run_button)
        self.recover_run_button = QPushButton("Recover Folder...")
        self.recover_run_button.setToolTip(
            "Rebuild a cautious result view from stage logs and artifacts when the JSON report is missing."
        )
        self.recover_run_button.clicked.connect(self._open_run_directory)
        buttons.addWidget(self.recover_run_button)
        buttons.addStretch(1)
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._start)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.stop_button)
        buttons.addWidget(close_button)
        content_layout.addLayout(buttons)

    def _set_subtree_state(self, item: QTreeWidgetItem, state: Qt.CheckState) -> None:
        item.setCheckState(0, state)
        for index in range(item.childCount()):
            self._set_subtree_state(item.child(index), state)

    def _refresh_ancestor_states(self, item: QTreeWidgetItem | None) -> None:
        while item is not None:
            states = [
                item.child(index).checkState(0)
                for index in range(item.childCount())
            ]
            if states and all(state == Qt.CheckState.Checked for state in states):
                item.setCheckState(0, Qt.CheckState.Checked)
            elif states and all(state == Qt.CheckState.Unchecked for state in states):
                item.setCheckState(0, Qt.CheckState.Unchecked)
            else:
                item.setCheckState(0, Qt.CheckState.PartiallyChecked)
            item = item.parent()

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        self.case_tree.blockSignals(True)
        for index in range(self.case_tree.topLevelItemCount()):
            self._set_subtree_state(self.case_tree.topLevelItem(index), state)
        self.case_tree.blockSignals(False)

    def _update_group_state(self, item: QTreeWidgetItem, column: int) -> None:
        if column != 0:
            return
        self.case_tree.blockSignals(True)
        if item.childCount():
            state = item.checkState(0)
            if state == Qt.CheckState.PartiallyChecked:
                state = Qt.CheckState.Unchecked
            for index in range(item.childCount()):
                self._set_subtree_state(item.child(index), state)
            item.setCheckState(0, state)
        self._refresh_ancestor_states(item.parent())
        self.case_tree.blockSignals(False)

    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        selected = (
            [(self._mesh_case_id, self._mesh_file.name)]
            if self._mesh_file is not None else [
                (case_id, item.text(0))
                for case_id, item in self._case_items.items()
                if item.checkState(0) == Qt.CheckState.Checked
            ]
        )
        if not selected:
            QMessageBox.information(self, "No mesh Cases selected", "Select at least one mesh Case.")
            return
        executable_text = self.executable_combo.currentText().strip()
        executable = Path(executable_text).expanduser().resolve() if executable_text else None
        requires_oneflow = any(
            mesh_case_requires_oneflow_executable(self._repository_root, case_id)
            for case_id, _ in selected
        )
        if requires_oneflow and (executable is None or not executable.is_file()):
            QMessageBox.warning(
                self, "OneFLOW executable required",
                "At least one selected Case includes OneFLOW stages. Select a valid executable.",
            )
            return
        if executable is not None and not executable.is_file():
            executable = None
        commands = []
        missing_mesh_inputs: list[str] = []
        for case_id, title in selected:
            arguments = [
                sys.executable, "-m", "tools.oneflow_gui.mesh_validation",
                "--mesh-case-id", case_id,
                "--cfd-template-id", str(self.template_combo.currentData()),
                "--max-steps", str(self.steps_spin.value()),
            ]
            if executable is not None:
                arguments.extend(["--oneflow-executable", str(executable)])
            if self._mesh_file is not None:
                if not self._mesh_file.is_file():
                    QMessageBox.warning(
                        self, "Mesh file required",
                        f"The current mesh no longer exists: {self._mesh_file}",
                    )
                    return
                arguments.extend(["--mesh-file", str(self._mesh_file)])
                if self._project_directory is not None:
                    arguments.extend(["--project-directory", str(self._project_directory)])
            else:
                existing_mesh_paths = existing_mesh_case_paths(
                    self._repository_root, case_id
                )
                if existing_mesh_paths is not None:
                    existing_mesh, project_directory = existing_mesh_paths
                    if not existing_mesh.is_file():
                        missing_mesh_inputs.append(f"{title}: {existing_mesh}")
                        item = self._case_items.get(case_id)
                        if item is not None:
                            item.setCheckState(0, Qt.CheckState.Unchecked)
                        continue
                    arguments.extend(["--mesh-file", str(existing_mesh)])
                    arguments.extend(["--project-directory", str(project_directory)])
            commands.append((f"Mesh-to-CFD: {title}", tuple(arguments), self._repository_root))
        if missing_mesh_inputs:
            QMessageBox.warning(
                self, "Generated mesh not available",
                "The following existing-mesh inputs were skipped. Generate them "
                "in Mesh Generation before selecting them:\n\n"
                + "\n".join(missing_mesh_inputs),
            )
        if not commands:
            return
        self.log_view.clear()
        self.result_tree.clear()
        self.artifact_tree.clear()
        self.raw_report_view.clear()
        self._reported_outcomes.clear()
        self.result_tabs.setCurrentIndex(0)
        self.result_summary.setText("Validation is running…")
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        log_directory = UserDataPaths.discover().mesh_validation / "_execution-logs"
        self._worker = MeshValidationWorker(commands, log_directory)
        self._worker.output_line.connect(self.log_view.appendPlainText)
        self._worker.report_ready.connect(self._report_ready)
        self._worker.completed.connect(self._completed)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()

    def _stop(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)

    def _open_run_report(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open Mesh Validation Report", "",
            "Mesh Validation Report (mesh-validation-report.json);;JSON files (*.json)",
        )
        if not filename:
            return
        report_path = Path(filename).expanduser().resolve()
        try:
            report = json.loads(report_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Cannot open validation report", str(error))
            return
        if not isinstance(report, dict) or report.get("schema") != "oneflow.run-report":
            QMessageBox.warning(
                self, "Unsupported validation report",
                "The selected file is not a supported oneflow.run-report JSON document.",
            )
            return
        artifacts = report.setdefault("artifacts", {})
        artifacts["report"] = str(report_path)
        artifacts["run_directory"] = str(report_path.parent)
        report["run_directory"] = str(report_path.parent)
        label = str(report.get("mesh_case_id") or report.get("workflow") or report_path.parent.name)
        self.result_tree.clear()
        self.artifact_tree.clear()
        self.log_view.clear()
        self._reported_outcomes.clear()
        self._report_ready(label, report)
        self.raw_report_view.setPlainText(json.dumps(report, indent=2, ensure_ascii=False))
        execution = report.get("execution") or {}
        console_value = execution.get("console_log") or (report.get("artifacts") or {}).get("console_log")
        if isinstance(console_value, str):
            console_path = Path(console_value).expanduser()
            if not console_path.is_absolute():
                console_path = report_path.parent / console_path
            if console_path.is_file():
                try:
                    self.log_view.setPlainText(console_path.read_text(encoding="utf-8", errors="replace"))
                except OSError as error:
                    self.log_view.setPlainText(f"Unable to read console log: {error}")
            else:
                self.log_view.setPlainText(f"Console log is missing: {console_path}")
        else:
            self.log_view.setPlainText("Console log: Not recorded in this report.")
        self.result_tabs.setCurrentIndex(0)

    def _open_run_directory(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Recover Mesh Validation Run", "")
        if directory:
            self._show_run_directory(Path(directory).expanduser().resolve())

    def _refresh_run_history(self) -> None:
        self.history_tree.clear()
        try:
            records = list_mesh_validation_runs(self._history_workspace)
        except (OSError, ValueError) as error:
            self.result_summary.setText(f"Unable to read run history: {error}")
            return
        for record in records:
            item = QTreeWidgetItem(self.history_tree, [
                record["mesh_case_id"], record["status"], record["started_at"], record["directory"],
            ])
            item.setData(0, Qt.ItemDataRole.UserRole, record["directory"])
            item.setToolTip(3, record["directory"])
        if not records:
            QTreeWidgetItem(self.history_tree, ["No validation runs", "", "", str(self._history_workspace)])

    def _open_selected_history_run(self, item: QTreeWidgetItem, column: int) -> None:
        directory = item.data(0, Qt.ItemDataRole.UserRole)
        if directory:
            self._show_run_directory(Path(directory))

    def _show_run_directory(self, run_directory: Path) -> None:
        report_path = run_directory / "mesh-validation-report.json"
        try:
            if report_path.is_file():
                report = json.loads(report_path.read_text(encoding="utf-8"))
                if not isinstance(report, dict) or report.get("schema") != "oneflow.run-report":
                    raise ValueError("The run report has an unsupported schema.")
                artifacts = report.setdefault("artifacts", {})
                artifacts["report"] = str(report_path)
                artifacts["run_directory"] = str(run_directory)
                report["run_directory"] = str(run_directory)
            else:
                report = recover_mesh_validation_run(run_directory)
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Cannot open validation run", str(error))
            return

        self.result_tree.clear()
        self.artifact_tree.clear()
        self.log_view.clear()
        self._reported_outcomes.clear()
        self._report_ready(run_directory.name, report)

        console_path = run_directory / "mesh-validation-console.log"
        if console_path.is_file():
            try:
                log_text = console_path.read_text(encoding="utf-8", errors="replace")
            except OSError as error:
                log_text = f"Unable to read recovered console log: {error}"
        else:
            sections = []
            for case_id in (
                "stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
                "oneflow_cfd_solve",
            ):
                stage_log = run_directory / case_id / "results" / "mesh-validation.log"
                if stage_log.is_file():
                    try:
                        sections.append(
                            f"===== {case_id} =====\n"
                            + stage_log.read_text(encoding="utf-8", errors="replace")
                        )
                    except OSError as error:
                        sections.append(f"===== {case_id} =====\nUnable to read log: {error}")
            log_text = (
                "\n\n".join(sections)
                if sections else "No durable console or stage logs were found."
            )
        self.log_view.setPlainText(log_text)
        self.result_tabs.setCurrentIndex(0)
        self._refresh_run_history()

    def _populate_artifacts(self, report: dict) -> None:
        artifacts = report.get("artifacts") or {}
        if not isinstance(artifacts, dict):
            artifacts = {}
        report_value = artifacts.get("report")
        run_directory = report.get("run_directory") or artifacts.get("run_directory")
        if isinstance(report_value, str):
            base_directory = Path(report_value).expanduser().parent
        elif isinstance(run_directory, str):
            base_directory = Path(run_directory).expanduser()
        else:
            base_directory = None

        def add_entry(parent, name: str, value) -> None:
            if isinstance(value, dict):
                item = QTreeWidgetItem(parent, [name, "Recorded", ""])
                for child_name, child_value in value.items():
                    add_entry(item, str(child_name), child_value)
                item.setExpanded(True)
                return
            if isinstance(value, list):
                item = QTreeWidgetItem(parent, [name, "Recorded", str(len(value)) + " entries"])
                for index, child_value in enumerate(value):
                    add_entry(item, "[" + str(index) + "]", child_value)
                item.setExpanded(True)
                return
            details = str(value)
            status = "Recorded"
            if isinstance(value, str) and (
                name.lower().endswith(("path", "file", "directory"))
                or "/" in value or "\\" in value
            ):
                candidate = Path(value).expanduser()
                if not candidate.is_absolute() and base_directory is not None:
                    candidate = base_directory / candidate
                try:
                    status = "Available" if candidate.exists() else "Missing"
                except OSError:
                    status = "Unresolved"
            item = QTreeWidgetItem(parent, [name, status, details])
            item.setToolTip(2, details)

        # Keep scalar run metadata visible, but render the file manifest as a
        # categorized inventory rather than a generic nested JSON list.
        for name, value in artifacts.items():
            if name != "files":
                add_entry(self.artifact_tree, str(name), value)

        manifest = artifacts.get("files")
        if isinstance(manifest, list):
            groups = {}
            for record in manifest:
                if not isinstance(record, dict) or not isinstance(record.get("path"), str):
                    continue
                category = str(record.get("category") or "Other")
                group = groups.get(category)
                if group is None:
                    group = QTreeWidgetItem(self.artifact_tree, [category, "", ""])
                    group.setExpanded(True)
                    groups[category] = group
                relative = Path(record["path"])
                status, candidate = _artifact_file_state(record["path"], base_directory)
                size = record.get("size_bytes")
                digest = record.get("sha256")
                details = record["path"]
                if isinstance(size, int):
                    details += " | " + format(size, ",") + " bytes"
                if isinstance(digest, str):
                    details += " | SHA-256: " + digest
                item = QTreeWidgetItem(
                    group, [str(record.get("name") or relative.name), status, details]
                )
                item.setToolTip(2, details)
                if candidate is not None:
                    item.setToolTip(0, str(candidate))
        elif manifest is None:
            # Backward compatibility: pre-inventory reports use legacy fields.
            pass
        elif manifest:
            add_entry(self.artifact_tree, "Unrecognized file manifest", manifest)

        lineage = report.get("geometry_lineage")
        if lineage:
            add_entry(self.artifact_tree, "Geometry lineage", lineage)
        if self.artifact_tree.topLevelItemCount() == 0:
            QTreeWidgetItem(
                self.artifact_tree,
                ["Artifacts", "Not recorded", "This report contains no artifact manifest."],
            )

    def _report_ready(self, label: str, report: dict) -> None:
        status = str(report.get("status", "unknown")).upper()
        exit_code = report.get("exit_code", "unknown")
        self.raw_report_view.setPlainText(json.dumps(report, indent=2, ensure_ascii=False))
        self._populate_artifacts(report)
        self._reported_outcomes.append(status)
        root = QTreeWidgetItem(self.result_tree, [label, status, f"Exit code: {exit_code}"])
        for stage_name, stage_status, stage_message in run_report_stage_rows(report):
            QTreeWidgetItem(root, [stage_name, stage_status, stage_message])
        lineage_item = QTreeWidgetItem(root, ["Geometry lineage", "", "Source and derived artifact trace"])
        for name, lineage_status, details in geometry_lineage_rows(report):
            QTreeWidgetItem(lineage_item, [name, lineage_status, details])
        lineage_item.setExpanded(True)
        root.setExpanded(True)
        failed = sum(item == "FAILED" for item in self._reported_outcomes)
        passed = sum(item == "SUCCESS" for item in self._reported_outcomes)
        self.result_summary.setText(
            f"Run reports: {len(self._reported_outcomes)}    "
            f"SUCCESS: {passed}    FAILED: {failed}"
        )

    def _completed(self, passed: bool, message: str) -> None:
        self.log_view.appendPlainText("\n" + message)
        verdict = "SUCCESS" if passed else "FAILED / CANCELLED"
        self.log_view.appendPlainText("OVERALL: " + verdict)
        self.result_summary.setText("OVERALL: " + verdict)
        self._refresh_run_history()

    def _worker_finished(self) -> None:
        self._worker = None
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._stop()
            event.ignore()
            return
        event.accept()
