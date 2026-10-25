"""Qt task dialog that launches the repository's existing test commands."""

from __future__ import annotations

import subprocess
import threading
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
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QAbstractItemView,
    QHeaderView,
    QVBoxLayout,
)

from .executable_discovery import executable_candidates
from .case_repository import load_cases
from .models import CaseDefinition
from .runner import create_run_id
from .user_data_paths import UserDataPaths
from .task_types import TaskRunResult, TaskRunStatus, TaskType, task_definition
from .test_runner import (
    TestCommand,
    available_suites,
    build_cmake_command,
    build_default_test_commands,
    build_regression_command,
    stage_regression_tests,
)


class TestCommandWorker(QThread):
    """Execute existing test commands sequentially and stream their output."""

    output_line = Signal(str)
    completed = Signal(bool, str)
    command_completed = Signal(object)

    def __init__(self, commands: list[TestCommand]) -> None:
        super().__init__()
        self._commands = commands
        self._cancel = threading.Event()
        self._process: subprocess.Popen[str] | None = None
        self._process_lock = threading.Lock()

    def request_stop(self) -> None:
        self._cancel.set()
        with self._process_lock:
            process = self._process
        if process is not None and process.poll() is None:
            process.terminate()

    def run(self) -> None:
        failed_labels: list[str] = []
        for command in self._commands:
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
            self.output_line.emit("$ " + subprocess.list2cmdline(list(command.arguments)))
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
                ) as process:
                    with self._process_lock:
                        self._process = process
                    if process.stdout is not None:
                        for line in process.stdout:
                            self.output_line.emit(line.rstrip("\r\n"))
                            if self._cancel.is_set() and process.poll() is None:
                                process.terminate()
                    return_code = process.wait()
                    with self._process_lock:
                        self._process = None
            except OSError as error:
                self.output_line.emit(f"ERROR: {error}")
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
                break
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
        if failed_labels:
            self.completed.emit(False, "Failed: " + ", ".join(failed_labels))
        else:
            self.completed.emit(True, f"Passed all {len(self._commands)} command(s).")


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
        self._database_cases: list[CaseDefinition] = load_cases(self._repository_root)
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

        unit_test_task = task_definition(TaskType.UNIT_TESTS)
        regression_task = task_definition(TaskType.CFD_GRID_REGRESSION)
        self.unit_tests_check = QCheckBox(
            f"{unit_test_task.title} (CTest and Python tests)"
        )
        self.unit_tests_check.setChecked(True)
        layout.addWidget(self.unit_tests_check)

        self.regression_check = QCheckBox(
            f"Optional {regression_task.title}"
        )
        self.regression_check.setChecked(False)
        layout.addWidget(self.regression_check)

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
        self.case_database_check = QCheckBox("Run Development Tests")
        self.case_database_check.setChecked(False)
        database_layout.insertWidget(0, self.case_database_check)
        self.case_database_button = QPushButton("Select cases...")
        self.case_database_button.setToolTip(
            "Choose which development cases to run. This does not start a run."
        )
        self.case_database_button.clicked.connect(self._open_case_database_tests)
        database_layout.addWidget(self.case_database_button)
        layout.addLayout(database_layout)

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

    def _open_case_database_tests(self) -> None:
        dialog = CaseDatabaseSelectionDialog(
            self._repository_root, self._database_cases, self
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._database_cases = dialog.selected_cases()
        count = len(self._database_cases)
        self.case_database_button.setText(
            f"Select cases... ({count} selected)" if count else "Select cases..."
        )

    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        if self._case_worker is not None and self._case_worker.isRunning():
            return
        commands: list[TestCommand] = []
        if self.unit_tests_check.isChecked():
            commands.extend(
                build_default_test_commands(
                    self._repository_root,
                    build_config=self._build_config,
                )
            )
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
                expected_build_executable = (
                    self._repository_root
                    / "build"
                    / "bin"
                    / self._build_config
                    / "OneFLOW.exe"
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
        run_database = self.case_database_check.isChecked()
        if run_database and not self._database_cases:
            QMessageBox.information(
                self, "No cases configured",
                "Use Select cases... to change the Development Tests selection.",
            )
            return
        if not commands and not run_database:
            QMessageBox.information(
                self, "No tests selected",
                "Select unit tests, CFD regression, or Development Tests.",
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

    def _launch_database_cases(self) -> None:
        executable = self.executable_combo.currentText().strip()
        if not executable or not Path(executable).expanduser().is_file():
            QMessageBox.warning(
                self, "Executable required",
                "Select a valid OneFLOW.exe before running Development Tests.",
            )
            return
        from .main_window import RunQueueWorker

        workspace_root = UserDataPaths.discover().grid_layout_test
        run_id = create_run_id()
        self.log_view.appendPlainText(
            f"\n=== Development Tests ===\n"
            f"Executable: {Path(executable).resolve()}\n"
            f"Workspace: {workspace_root}"
        )
        self.run_button.setEnabled(False)
        self.build_button.setEnabled(False)
        self.rebuild_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._case_worker = RunQueueWorker(
            self._database_cases,
            Path(executable).expanduser().resolve(),
            workspace_root,
            run_id,
            execution_mode="case_database",
        )
        self._case_worker.log_line.connect(self.log_view.appendPlainText)
        self._case_worker.case_result.connect(self._database_case_result)
        self._case_worker.finished.connect(self._database_worker_finished)
        self._case_worker.start()

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
        self._worker = TestCommandWorker(commands)
        self._worker.output_line.connect(self.log_view.appendPlainText)
        self._worker.completed.connect(self._completed)
        self._worker.command_completed.connect(self._command_completed)
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

    def _completed(self, passed: bool, message: str) -> None:
        self.log_view.appendPlainText(f"\n{message}")
        if self._run_database_after_commands and passed:
            self.log_view.appendPlainText("Unit/regression commands passed; starting Development Tests.")
            return
        self.log_view.appendPlainText("FINAL STATUS: PASS" if passed else "FINAL STATUS: FAIL")

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
        repository_root: Path,
        selected_cases: list[CaseDefinition],
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._repository_root = repository_root
        self._cases = load_cases(repository_root)
        selected_ids = {case.case_id for case in selected_cases}
        self.setWindowTitle("Select Development Test Cases")
        self.resize(760, 560)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Select development cases to include when you press Run. "
            "This dialog only saves the selection; it does not start calculations."
        ))
        self.case_table = QTableWidget(0, 3)
        self.case_table.setHorizontalHeaderLabels(["Run", "Case", "Configuration"])
        self.case_table.verticalHeader().setVisible(False)
        self.case_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.case_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.case_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        for case in self._cases:
            row = self.case_table.rowCount()
            self.case_table.insertRow(row)
            check = QTableWidgetItem()
            check.setFlags(check.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            check.setCheckState(
                Qt.CheckState.Checked if case.case_id in selected_ids else Qt.CheckState.Unchecked
            )
            self.case_table.setItem(row, 0, check)
            self.case_table.setItem(row, 1, QTableWidgetItem(case.title))
            self.case_table.setItem(
                row, 2, QTableWidgetItem("Typed profile" if case.profile else "Legacy config")
            )
        layout.addWidget(self.case_table, 1)
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

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for row in range(self.case_table.rowCount()):
            self.case_table.item(row, 0).setCheckState(state)

    def selected_cases(self) -> list[CaseDefinition]:
        selected_ids = {
            case.case_id
            for row, case in enumerate(self._cases)
            if self.case_table.item(row, 0).checkState() == Qt.CheckState.Checked
        }
        return [case for case in self._cases if case.case_id in selected_ids]
