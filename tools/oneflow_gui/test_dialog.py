"""Qt task dialog that launches the repository's existing test commands."""

from __future__ import annotations

import shutil
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QSpinBox,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QHeaderView,
    QVBoxLayout,
)

from .executable_discovery import executable_candidates
from .case_repository import (
    load_cases,
    select_cases_by_execution_kind,
    select_cases_for_execution,
    select_cases_for_selection_dialog,
)
from .models import CaseDefinition
from .process_control import process_group_options, terminate_and_reap, terminate_process_tree
from .runner import prepare_case_database_workspace
from .user_data_paths import UserDataPaths
from .task_types import TaskRunResult, TaskRunStatus, TaskType, task_definition
from .test_runner import (
    TestCommand,
    available_suites,
    build_cmake_command,
    build_ctest_test_commands,
    build_python_test_commands,
    build_regression_command,
    format_command_arguments,
    format_final_test_status,
    oneflow_executable_path,
    stage_regression_tests,
    summarize_test_output,
)
from .test_diagnostics import format_failure_summary, parse_test_diagnostics
from .execution_capture import (
    capture_python_unittest_execution,
    is_python_unittest_command,
)


class TestCommandWorker(QThread):
    """Execute existing test commands sequentially and stream their output."""

    output_line = Signal(str)
    completed = Signal(bool, str)
    command_completed = Signal(object)
    diagnostics = Signal(object)

    def __init__(self, commands: list[TestCommand], *, build_config: str | None = None) -> None:
        super().__init__()
        self._commands = commands
        self._build_config = build_config
        self._cancel = threading.Event()
        self._process: subprocess.Popen[str] | None = None
        self._process_lock = threading.Lock()

    def request_stop(self) -> None:
        self._cancel.set()
        with self._process_lock:
            process = self._process
        if process is not None:
            terminate_process_tree(process)

    def run(self) -> None:
        failed_labels: list[str] = []
        passed_labels: list[str] = []
        run_output: list[str] = []
        failed_diagnostics = []
        for command in self._commands:
            command_output: list[str] = []
            if self._cancel.is_set():
                self.command_completed.emit(
                    TaskRunResult(
                        task_type=command.task_type,
                        status=TaskRunStatus.CANCELLED,
                        return_code=None,
                        execution_directory=command.working_directory,
                        message="Test run cancelled before command start.",
                    )
                )
                self.completed.emit(False, "Test run cancelled.")
                return
            self.output_line.emit(f"\n=== {command.label} ===")
            self.output_line.emit("$ " + format_command_arguments(command.arguments))
            started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
            try:
                if command.label.startswith("CMake configure"):
                    command.working_directory.mkdir(parents=True, exist_ok=True)
                with subprocess.Popen(
                    list(command.arguments),
                    cwd=command.working_directory,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    bufsize=1,
                    **process_group_options(),
                ) as process:
                    with self._process_lock:
                        self._process = process
                    if self._cancel.is_set():
                        terminate_and_reap(process)
                    try:
                        if process.stdout is not None:
                            for line in process.stdout:
                                output_line = line.rstrip("\r\n")
                                command_output.append(output_line)
                                self.output_line.emit(output_line)
                                if self._cancel.is_set():
                                    terminate_and_reap(process)
                                    break
                        return_code = process.wait()
                    finally:
                        with self._process_lock:
                            if self._process is process:
                                self._process = None
            except OSError as error:
                self.output_line.emit(f"ERROR: {error}")
                run_output.append(f"ERROR: {error}")
                self.command_completed.emit(
                    TaskRunResult(
                        task_type=command.task_type,
                        status=TaskRunStatus.FAILED,
                        return_code=None,
                        execution_directory=command.working_directory,
                        message=str(error),
                    )
                )
                failed_labels.append(command.label)
                run_output.extend(command_output)
                break
            finished_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
            if is_python_unittest_command(command):
                try:
                    user_paths = UserDataPaths.discover()
                    record = capture_python_unittest_execution(
                        command,
                        repository_root=command.working_directory,
                        records_path=user_paths.verification_execution_records,
                        artifact_directory=user_paths.test_run_logs,
                        started_at=started_at,
                        finished_at=finished_at,
                        exit_code=return_code,
                        output_lines=command_output,
                        cancelled=self._cancel.is_set(),
                        configuration=self._build_config or "",
                    )
                    self.output_line.emit(
                        f"EXECUTION RECORD: {record.record_id} / {record.outcome}"
                    )
                except Exception as error:
                    self.output_line.emit(
                        f"WARNING: unable to persist execution evidence: {error}"
                    )
            run_output.extend(command_output)
            if self._cancel.is_set():
                self.command_completed.emit(
                    TaskRunResult(
                        task_type=command.task_type,
                        status=TaskRunStatus.CANCELLED,
                        return_code=return_code,
                        execution_directory=command.working_directory,
                        message="Test run cancelled.",
                    )
                )
                self.completed.emit(False, "Test run cancelled.")
                return
            if return_code != 0:
                configuration = self._build_config
                for index, argument in enumerate(command.arguments[:-1]):
                    if argument == "-C":
                        configuration = command.arguments[index + 1]
                        break
                    if argument.startswith("CMAKE_BUILD_TYPE="):
                        configuration = argument.split("=", 1)[1]
                framework = "GoogleTest" if "CTest" in command.label else (
                    "unittest" if "unittest" in command.arguments else ""
                )
                if framework:
                    records = parse_test_diagnostics(
                        command_output,
                        task_label=command.label,
                        framework=framework,
                        configuration=configuration,
                        working_directory=command.working_directory,
                        started_at=started_at,
                        finished_at=finished_at,
                        exit_code=return_code,
                    )
                    if records:
                        failed_diagnostics.extend(records)
                        self.diagnostics.emit(records)
                self.command_completed.emit(
                    TaskRunResult(
                        task_type=command.task_type,
                        status=TaskRunStatus.FAILED,
                        return_code=return_code,
                        execution_directory=command.working_directory,
                        message=f"{command.label} failed with exit code {return_code}.",
                    )
                )
                failed_labels.append(command.label)
                self.output_line.emit(f"FAILED: exit code {return_code}")
                break
            passed_labels.append(command.label)
            self.command_completed.emit(
                TaskRunResult(
                    task_type=command.task_type,
                    status=TaskRunStatus.SUCCEEDED,
                    return_code=return_code,
                    execution_directory=command.working_directory,
                    message=f"{command.label} passed.",
                )
            )
            self.output_line.emit("PASSED")

        not_run = len(self._commands) - len(passed_labels) - len(failed_labels)
        report = [
            "Test command summary:",
            f"  Total: {len(self._commands)}",
            f"  Passed: {len(passed_labels)}",
            f"  Failed: {len(failed_labels)}",
            f"  Not run: {not_run}",
        ]
        report.extend(
            summarize_test_output(
                run_output, include_failed_cases=not failed_diagnostics
            )
        )
        report.extend(format_failure_summary(failed_diagnostics))
        if failed_labels:
            report.insert(1, "  Failed commands: " + ", ".join(failed_labels))
            self.completed.emit(False, "\n".join(report))
        else:
            self.completed.emit(True, "\n".join(report))


class TestRunnerDialog(QDialog):
    """Run unit tests, CFD regression, or catalog cases from one test entry point."""

    def __init__(
        self,
        repository_root: Path,
        parent=None,
        *,
        build_config: str = "Debug",
        executable: str | None = None,
    ) -> None:
        super().__init__(parent)
        if build_config not in {"Debug", "Release"}:
            raise ValueError(f"Unsupported build configuration: {build_config!r}")
        self._repository_root = Path(repository_root).resolve()
        self._build_config = build_config
        self._worker: TestCommandWorker | None = None
        self._case_worker = None
        all_database_cases = load_cases(self._repository_root)
        self._all_database_cases = all_database_cases
        self._available_database_cases_by_kind: dict[str, list[CaseDefinition]] = {
            "oneflow": select_cases_by_execution_kind(all_database_cases, "oneflow"),
            "python": select_cases_by_execution_kind(all_database_cases, "python"),
        }
        self._database_cases_by_kind = {
            kind: list(cases)
            for kind, cases in self._available_database_cases_by_kind.items()
        }
        self._run_database_after_commands = False
        self.setWindowTitle("OneFLOW Tests")
        self.resize(900, 650)

        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "Run unit tests, CFD/grid regression, or catalog case tests. "
                "Development Tests run selected catalog cases from a refreshed, isolated test workspace."
            )
        )

        regression_task = task_definition(TaskType.CFD_GRID_REGRESSION)
        self.ctest_check = QCheckBox("C++ Unit Tests (CTest)")
        self.ctest_check.setToolTip(
            "Run C++ GoogleTest targets registered in the configured CMake build tree."
        )
        self.ctest_check.setChecked(True)
        layout.addWidget(self.ctest_check)
        self.python_tests_check = QCheckBox("Python Unit Tests (unittest)")
        self.python_tests_check.setToolTip(
            "Run the explicitly registered Python unittest modules with the selected Python interpreter."
        )
        self.python_tests_check.setChecked(True)
        layout.addWidget(self.python_tests_check)

        self.regression_check = QCheckBox(
            f"Optional {regression_task.title}"
        )
        self.regression_check.setChecked(False)
        layout.addWidget(self.regression_check)

        self.gmsh_case_check = QCheckBox(
            "Experimental Gmsh-to-CGNS Case (manual)"
        )
        self.gmsh_case_check.setChecked(False)
        self.gmsh_case_check.setToolTip(
            "Run the experimental Gmsh MSH 4.1 Case (or legacy 2.x) from a temporary copy. "
            "The generated .cgns file is an explicitly marked placeholder, not a CGNS file."
        )
        layout.addWidget(self.gmsh_case_check)

        mesh_validation_layout = QHBoxLayout()
        mesh_validation_description = QLabel(
            "Run complete STL-to-CFD validation in its own workspace."
        )
        mesh_validation_description.setWordWrap(True)
        mesh_validation_layout.addWidget(mesh_validation_description, 1)
        self.mesh_validation_button = QPushButton("Mesh-to-CFD Validation...")
        self.mesh_validation_button.setToolTip(
            "Open the dedicated validation workbench to select multiple mesh Cases, "
            "a CFD consumer template, and solver smoke-test settings."
        )
        self.mesh_validation_button.clicked.connect(self._open_mesh_validation)
        mesh_validation_layout.addWidget(self.mesh_validation_button)
        layout.addLayout(mesh_validation_layout)

        form = QFormLayout()
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(42)
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if executable:
            self.executable_combo.setEditText(executable)
        elif self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        self.suite_combo = QComboBox()
        self.suite_combo.addItems(available_suites(self._repository_root))
        default_index = self.suite_combo.findText("ci-serial.txt")
        if default_index >= 0:
            self.suite_combo.setCurrentIndex(default_index)
        form.addRow("OneFLOW executable", self.executable_combo)
        form.addRow("Existing suite", self.suite_combo)
        layout.addLayout(form)

        database_layout = QHBoxLayout()
        database_description = QLabel(
            "Run all selected development cases from copies staged under the fixed test workspace."
        )
        database_description.setWordWrap(True)
        database_layout.addWidget(database_description, 1)
        self.case_database_check = QCheckBox("Run Development Cases")
        self.case_database_check.setChecked(False)
        database_layout.insertWidget(0, self.case_database_check)
        self.python_development_check = QCheckBox("Python reference cases")
        self.python_development_check.setChecked(False)
        self.python_development_check.setToolTip(
            "Runs Python reference cases. Together with Run Development Cases, "
            "you can run Python only, native C++ only, or both. Selections are retained by type."
        )
        self.python_development_check.toggled.connect(
            self._update_case_database_button_text
        )
        self.case_database_check.toggled.connect(
            self._update_case_database_button_text
        )
        database_layout.insertWidget(1, self.python_development_check)
        self.case_database_button = QPushButton("Select cases...")
        self.case_database_button.setToolTip(
            "Choose which development cases to run. This does not start a run."
        )
        self.case_database_button.clicked.connect(self._open_case_database_tests)
        database_layout.addWidget(self.case_database_button)
        layout.addLayout(database_layout)
        self._update_case_database_button_text()

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.log_view.setPlaceholderText("Test output will appear here.")
        layout.addWidget(self.log_view, 1)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        self.build_button = QPushButton("Build")
        self.build_button.setToolTip("Build the selected CMake configuration.")
        self.build_button.clicked.connect(lambda: self._start_build(rebuild=False))
        self.rebuild_button = QPushButton("Rebuild")
        self.rebuild_button.setToolTip("Clean and rebuild the selected CMake configuration.")
        self.rebuild_button.clicked.connect(lambda: self._start_build(rebuild=True))
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._start)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(self.build_button)
        buttons.addWidget(self.rebuild_button)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.stop_button)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)

    def _open_mesh_validation(self) -> None:
        from .mesh_validation_dialog import MeshValidationDialog

        dialog = MeshValidationDialog(
            self._repository_root,
            executable=self.executable_combo.currentText().strip(),
            parent=self,
        )
        dialog.exec()

    def _development_case_selection_kinds(self) -> tuple[str, ...]:
        if self.case_database_check.isChecked() and self.python_development_check.isChecked():
            return ("oneflow", "python")
        if self.python_development_check.isChecked():
            return ("python",)
        return ("oneflow",)

    def _update_case_database_button_text(self, *_args) -> None:
        kinds = self._development_case_selection_kinds()
        labels = []
        for kind in kinds:
            label = "Python" if kind == "python" else "C++"
            labels.append(f"{label} ({len(self._database_cases_by_kind[kind])})")
        self.case_database_button.setText("Select cases... [" + " + ".join(labels) + "]")

    def _open_case_database_tests(self) -> None:
        available_cases = select_cases_for_selection_dialog(
            self._all_database_cases,
            run_native=self.case_database_check.isChecked(),
            run_python=self.python_development_check.isChecked(),
        )
        kinds = {case.execution_kind for case in available_cases}
        selected_cases = [
            case for kind in kinds for case in self._database_cases_by_kind[kind]
        ]
        dialog = CaseDatabaseSelectionDialog(available_cases, selected_cases, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        selected = dialog.selected_cases()
        for kind in kinds:
            self._database_cases_by_kind[kind] = [
                case for case in selected if case.execution_kind == kind
            ]
        self._update_case_database_button_text()

    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        if self._case_worker is not None and self._case_worker.isRunning():
            return
        commands: list[TestCommand] = []
        if self.ctest_check.isChecked():
            ctest_commands = build_ctest_test_commands(
                self._repository_root,
                build_config=self._build_config,
            )
            if not ctest_commands:
                QMessageBox.warning(
                    self,
                    "CTest is not configured",
                    "No build/CTestTestfile.cmake was found. Configure the CMake build tree first.",
                )
                return
            commands.extend(ctest_commands)
        if self.python_tests_check.isChecked():
            commands.extend(build_python_test_commands(self._repository_root))
        if self.regression_check.isChecked():
            executable = self.executable_combo.currentText().strip()
            suite_name = self.suite_combo.currentText().strip()
            try:
                if not suite_name:
                    raise ValueError("No regression suite is available.")
                staged_test_directory = stage_regression_tests(
                    self._repository_root,
                    UserDataPaths.discover().cfd_regression_workspace,
                )
                executable_path = Path(executable).expanduser().resolve()
                expected_build_executable = oneflow_executable_path(
                    self._repository_root,
                    self._build_config,
                ).resolve()
                executable_will_be_built = (
                    executable_path == expected_build_executable
                    and any(command.label.startswith("CMake build") for command in commands)
                )
                commands.append(
                    build_regression_command(
                        self._repository_root,
                        executable,
                        suite_name,
                        test_directory=staged_test_directory,
                        validate_executable=not executable_will_be_built,
                    )
                )
            except (OSError, ValueError) as error:
                QMessageBox.warning(self, "Cannot start regression", str(error))
                return
        if self.gmsh_case_check.isChecked():
            executable = self.executable_combo.currentText().strip()
            try:
                if not executable:
                    raise ValueError("Select a OneFLOW executable for the Gmsh Case.")
                executable_path = Path(executable).expanduser().resolve()
                if not executable_path.is_file():
                    raise FileNotFoundError(
                        f"OneFLOW executable not found: {executable_path}"
                    )
                commands.append(
                    TestCommand(
                        "Experimental Gmsh-to-CGNS Case",
                        (
                            sys.executable,
                            "-m",
                            "tools.oneflow_gui.gmsh_case_runner",
                            "--repository-root",
                            str(self._repository_root),
                            "--executable",
                            str(executable_path),
                        ),
                        self._repository_root,
                        task_type=TaskType.CFD_GRID_REGRESSION,
                    )
                )
            except (OSError, ValueError) as error:
                QMessageBox.warning(self, "Cannot start Gmsh Case", str(error))
                return

        run_database = (
            self.case_database_check.isChecked()
            or self.python_development_check.isChecked()
        )
        if run_database and not self._selected_development_cases():
            mode = "native C++ OneFLOW and Python" if (
                self.case_database_check.isChecked() and self.python_development_check.isChecked()
            ) else ("Python" if self.python_development_check.isChecked() else "native C++ OneFLOW")
            QMessageBox.information(
                self, "No cases configured",
                f"No {mode} cases are selected. Use Select cases... to change the selection.",
            )
            return
        if not commands and not run_database:
            QMessageBox.information(
                self, "No tests selected",
                "Select CTest, Python Unit Tests, CFD regression, the experimental Gmsh Case, or Development Tests. Mesh-to-CFD Validation runs in its own workbench.",
            )
            return

        self._run_database_after_commands = bool(commands and run_database)
        if commands:
            self._launch_commands(commands)
        elif run_database:
            self._launch_database_cases()

    def _start_build(self, *, rebuild: bool) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        command = build_cmake_command(
            self._repository_root,
            self._build_config,
            rebuild=rebuild,
        )
        self._launch_commands([command])

    def _selected_development_cases(self) -> list[CaseDefinition]:
        selected = [
            case
            for cases in self._database_cases_by_kind.values()
            for case in cases
        ]
        return select_cases_for_execution(
            selected,
            run_native=self.case_database_check.isChecked(),
            run_python=self.python_development_check.isChecked(),
        )

    def _launch_database_cases(self) -> None:
        from .runner import build_command

        selected_cases = self._selected_development_cases()
        workspace_root = UserDataPaths.discover().grid_layout_test
        executable_text = self.executable_combo.currentText().strip()
        needs_oneflow = any(case.execution_kind == "oneflow" for case in selected_cases)
        executable = Path(executable_text).expanduser().resolve() if executable_text else None
        if needs_oneflow and (executable is None or not executable.is_file()):
            QMessageBox.warning(
                self, "Executable required",
                "The selected development cases include OneFLOW cases. Select a valid OneFLOW executable.",
            )
            return
        try:
            staged = prepare_case_database_workspace(selected_cases, workspace_root)
            commands: list[TestCommand] = []
            for case in selected_cases:
                case_directory = staged[case.case_id]
                if case.execution_kind == "python":
                    if case.execution_script is None:
                        raise ValueError(f"{case.case_id}: Python execution script is missing.")
                    arguments = (sys.executable, str(case.execution_script), "--case-directory", str(case_directory))
                    working_directory = self._repository_root
                else:
                    assert executable is not None
                    arguments = tuple(build_command(executable, case_directory))
                    working_directory = executable.parent
                commands.append(TestCommand(
                    f"Development case: {case.title}", arguments, working_directory,
                    task_type=TaskType.CASE_DATABASE,
                ))
        except (OSError, ValueError, shutil.Error) as error:
            QMessageBox.warning(self, "Cannot stage development cases", str(error))
            return
        if not commands:
            QMessageBox.information(self, "No cases selected", "Select at least one development case.")
            return
        enabled_modes = []
        if any(case.execution_kind == "oneflow" for case in selected_cases):
            enabled_modes.append("native C++ OneFLOW")
        if any(case.execution_kind == "python" for case in selected_cases):
            enabled_modes.append("Python reference")
        self.log_view.appendPlainText(
            f"Development case modes: {' + '.join(enabled_modes)}\n"
            f"Development case workspace: {workspace_root}\nStaged {len(commands)} selected case(s)."
        )
        self._launch_commands(commands)

    def _database_case_result(self, result) -> None:
        self.log_view.appendPlainText(
            f"[{result.case_id}] {result.state.value}: {result.message}\n"
            f"Workspace: {result.workspace}"
        )

    def _database_worker_finished(self) -> None:
        self._case_worker = None
        self.run_button.setEnabled(True)
        self.build_button.setEnabled(True)
        self.rebuild_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.log_view.appendPlainText("Development Tests finished.")

    def _launch_commands(self, commands: list[TestCommand]) -> None:
        self.log_view.clear()
        self.run_button.setEnabled(False)
        self.build_button.setEnabled(False)
        self.rebuild_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = TestCommandWorker(commands, build_config=self._build_config)
        self._worker.output_line.connect(self.log_view.appendPlainText)
        self._worker.completed.connect(self._completed)
        self._worker.command_completed.connect(self._command_completed)
        self._worker.diagnostics.connect(self._show_diagnostics)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()

    def _stop(self) -> None:
        stopping = False
        if self._worker is not None:
            self._worker.request_stop()
            stopping = True
        if self._case_worker is not None:
            self._case_worker.request_stop()
            stopping = True
        if stopping:
            self.stop_button.setEnabled(False)

    def _command_completed(self, result: TaskRunResult) -> None:
        self.log_view.appendPlainText(
            f"RESULT: {result.task_type.value} / {result.status.value} "
            f"/ exit={result.return_code} / cwd={result.execution_directory}"
        )

    def _show_diagnostics(self, diagnostics) -> None:
        for item in diagnostics:
            location = item.source_file or "unknown"
            if item.source_line is not None:
                location += f":{item.source_line}"
            cause = item.exception_type or "unknown"
            if item.message:
                cause += f": {item.message}"
            self.log_view.appendPlainText(
                "\nSTRUCTURED FAILURE"
                f"\n  Task: {item.task_label}"
                f"\n  Framework: {item.framework}"
                f"\n  Test: {item.test_id or 'unknown'}"
                f"\n  Location: {location}"
                f"\n  Cause: {cause}"
                f"\n  Configuration: {item.configuration or 'unknown'}"
                f"\n  Working directory: {item.working_directory or 'unknown'}"
                f"\n  Started: {item.started_at or 'unknown'}"
                f"\n  Finished: {item.finished_at or 'unknown'}"
                f"\n  Exit code: {item.exit_code}"
                f"\n  Contract links: {', '.join(item.contract_ids) or 'not resolved'}"
                f"\n  Verification impact: {item.verification_impact}"
            )

    def _completed(self, passed: bool, message: str) -> None:
        self.log_view.appendPlainText(f"\n{message}")
        if self._run_database_after_commands and passed:
            self.log_view.appendPlainText("Unit/regression commands passed; starting Development Tests.")
            return
        self.log_view.appendPlainText(format_final_test_status(passed, message))

    def _worker_finished(self) -> None:
        self.run_button.setEnabled(True)
        self.build_button.setEnabled(True)
        self.rebuild_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self._worker = None
        if self._run_database_after_commands:
            self._run_database_after_commands = False
            self._launch_database_cases()

    def closeEvent(self, event) -> None:
        command_running = self._worker is not None and self._worker.isRunning()
        cases_running = self._case_worker is not None and self._case_worker.isRunning()
        if command_running or cases_running:
            self._stop()
            event.ignore()
            return
        event.accept()


class CaseDatabaseSelectionDialog(QDialog):
    """Configure which development cases the parent Tests page should run."""

    def __init__(
        self,
        available_cases: list[CaseDefinition],
        selected_cases: list[CaseDefinition],
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._cases = list(available_cases)
        selected_ids = {case.case_id for case in selected_cases}
        self.setWindowTitle("Select Development Test Cases")
        self.resize(760, 560)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Select development cases to include when you press Run. "
            "This dialog only saves the selection; it does not start calculations."
        ))
        self.case_tree = QTreeWidget()
        self.case_tree.setColumnCount(3)
        self.case_tree.setHeaderLabels(["Run", "Case / Group", "Configuration"])
        self.case_tree.setRootIsDecorated(True)
        self.case_tree.setAlternatingRowColors(True)
        self.case_tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        dimension_items: dict[str, QTreeWidgetItem] = {}
        category_items: dict[tuple[str, str], QTreeWidgetItem] = {}
        for case in self._cases:
            dimension = f"{case.dimension}D" if case.dimension in (2, 3) else "Unclassified"
            dimension_item = dimension_items.get(dimension)
            if dimension_item is None:
                dimension_item = QTreeWidgetItem(self.case_tree)
                dimension_item.setText(1, dimension)
                dimension_item.setFlags(
                    dimension_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                )
                dimension_items[dimension] = dimension_item

            category_key = (dimension, case.group or "General")
            category_item = category_items.get(category_key)
            if category_item is None:
                category_item = QTreeWidgetItem(dimension_item)
                category_item.setText(1, category_key[1])
                category_item.setFlags(
                    category_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                )
                category_items[category_key] = category_item

            case_item = QTreeWidgetItem(category_item)
            case_item.setFlags(case_item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            case_item.setData(0, Qt.ItemDataRole.UserRole, case.case_id)
            case_item.setText(1, case.title)
            configuration = (
                f"Python: {case.execution_script.name}"
                if case.execution_kind == "python" and case.execution_script is not None
                else ("Typed profile" if case.profile else "Legacy config")
            )
            case_item.setText(2, configuration)
            case_item.setCheckState(
                0,
                Qt.CheckState.Checked
                if case.case_id in selected_ids
                else Qt.CheckState.Unchecked,
            )

        for category_item in category_items.values():
            self._refresh_ancestor_states(category_item)
        self.case_tree.itemChanged.connect(self._on_case_item_changed)
        self.case_tree.expandToDepth(1)
        layout.addWidget(self.case_tree, 1)
        actions = QHBoxLayout()
        select_all = QPushButton("Select all")
        select_all.clicked.connect(lambda: self._set_all(True))
        select_none = QPushButton("Select none")
        select_none.clicked.connect(lambda: self._set_all(False))
        actions.addWidget(select_all)
        actions.addWidget(select_none)
        actions.addStretch(1)
        ok_button = QPushButton("OK")
        ok_button.clicked.connect(self.accept)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        actions.addWidget(ok_button)
        actions.addWidget(cancel_button)
        layout.addLayout(actions)

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

    def _on_case_item_changed(self, item: QTreeWidgetItem, column: int) -> None:
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

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        self.case_tree.blockSignals(True)
        for index in range(self.case_tree.topLevelItemCount()):
            self._set_subtree_state(self.case_tree.topLevelItem(index), state)
        self.case_tree.blockSignals(False)

    def selected_cases(self) -> list[CaseDefinition]:
        selected_ids = set()
        def visit(item: QTreeWidgetItem) -> None:
            if item.childCount() == 0:
                case_id = item.data(0, Qt.ItemDataRole.UserRole)
                if case_id and item.checkState(0) == Qt.CheckState.Checked:
                    selected_ids.add(case_id)
                return
            for index in range(item.childCount()):
                visit(item.child(index))
        for index in range(self.case_tree.topLevelItemCount()):
            visit(self.case_tree.topLevelItem(index))
        return [case for case in self._cases if case.case_id in selected_ids]
