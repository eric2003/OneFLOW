"""Qt desktop interface for selecting and running OneFLOW cases."""

from __future__ import annotations

import threading
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .case_repository import load_cases
from .executable_discovery import executable_candidates
from .models import CaseDefinition, CaseRunResult, RunState
from .runner import create_run_id, run_case
from .user_data_paths import UserDataPaths
from .workflow_dialog import WorkflowDialog
from .test_dialog import TestCommandWorker, TestRunnerDialog
from .test_runner import build_cmake_command
from .run_history_dialog import RunHistoryDialog
from .workflow_repository import WorkflowRepository, WorkflowRunRecord
from .mesh_layout_dialog import MeshLayoutDialog


class RunQueueWorker(QThread):
    """Execute selected cases serially away from the UI thread."""

    log_line = Signal(str)
    case_state = Signal(str, str)
    case_result = Signal(object)

    def __init__(
        self,
        cases: list[CaseDefinition],
        executable: Path,
        workspace_root: Path,
        run_id: str,
        execution_mode: str = "default_workspace",
    ) -> None:
        super().__init__()
        self._cases = cases
        self._executable = executable
        self._workspace_root = workspace_root
        self._run_id = run_id
        self._execution_mode = execution_mode
        self._cancel_event = threading.Event()

    def request_stop(self) -> None:
        """Ask the active process to stop and prevent subsequent launches."""
        self._cancel_event.set()

    def run(self) -> None:
        for index, case in enumerate(self._cases):
            if self._cancel_event.is_set():
                result = CaseRunResult(
                    case.case_id,
                    RunState.CANCELLED,
                    None,
                    self._workspace_root,
                    "Cancelled before launch.",
                )
                self.case_result.emit(result)
                self.case_state.emit(case.case_id, RunState.CANCELLED.value)
                continue

            self.case_state.emit(case.case_id, RunState.RUNNING.value)
            self.log_line.emit(
                f"\n=== Case {index + 1}/{len(self._cases)}: {case.title} ==="
            )
            result = run_case(
                case,
                self._executable,
                self._workspace_root,
                self._run_id,
                self._cancel_event,
                self.log_line.emit,
                execution_mode=self._execution_mode,
            )
            self.case_result.emit(result)
            self.case_state.emit(case.case_id, result.state.value)


class MainWindow(QMainWindow):
    """Desktop shell; case loading and process execution live in services."""

    def __init__(
        self,
        repository_root: Path | None = None,
        *,
        single_case_mode: bool = False,
    ) -> None:
        super().__init__()
        self._single_case_mode = single_case_mode
        self._repository_root = (
            repository_root or Path(__file__).resolve().parents[2]
        ).resolve()
        self._cases: list[CaseDefinition] = []
        self._worker: RunQueueWorker | None = None
        self._task_worker: TestCommandWorker | None = None
        self._row_by_id: dict[str, int] = {}
        self._history_repository = WorkflowRepository(UserDataPaths.discover().workflow_data)
        self._active_case_run_record: WorkflowRunRecord | None = None
        self._active_case_run_steps: list[dict] = []
        self._active_case_run_log_path: Path | None = None

        self.setWindowTitle(
            "Case Database" if self._single_case_mode else "Case Workbench"
        )
        self.resize(1100, 760)
        self._build_ui()
        self._load_cases()

    def _build_ui(self) -> None:
        """Build a compact layout with configuration, queue, and log areas."""
        central = QWidget(self)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(20, 18, 20, 18)
        root_layout.setSpacing(12)

        heading = QLabel(
            "Case Database" if self._single_case_mode else "Case Workbench"
        )
        heading.setObjectName("pageTitle")
        subtitle = QLabel(
            "Select database cases and run them from protected workspace copies."
            if self._single_case_mode
            else "Select cases, run them in isolated workspaces, and review process logs."
        )
        subtitle.setObjectName("subtitle")
        root_layout.addWidget(heading)
        root_layout.addWidget(subtitle)

        executable_layout = QHBoxLayout()
        executable_layout.addWidget(QLabel("OneFLOW executable"))
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.executable_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.executable_combo.setMinimumContentsLength(38)
        self.executable_combo.setToolTip(
            "Detected paths are ordered by preference: repository Debug, repository "
            "Release, C:/dev, then conventional Windows installation folders. "
            "You can also type or browse to another executable."
        )
        self.executable_combo.lineEdit().setPlaceholderText(
            "Select or enter the path to OneFLOW.exe"
        )
        self._refresh_executable_candidates()
        executable_layout.addWidget(self.executable_combo, 1)
        executable_browse = QPushButton("Browse...")
        executable_browse.clicked.connect(self._choose_executable)
        executable_layout.addWidget(executable_browse)
        root_layout.addLayout(executable_layout)

        location_layout = QHBoxLayout()
        location_layout.addWidget(QLabel("Execution location"))
        self.execution_location_combo = QComboBox()
        self.execution_location_combo.addItem(
            "Copy to OneFLOW default run workspace", "default_workspace"
        )
        self.execution_location_combo.addItem(
            "Copy to specified run directory", "custom_workspace"
        )
        self.execution_location_combo.currentIndexChanged.connect(
            self._execution_location_changed
        )
        location_layout.addWidget(self.execution_location_combo, 1)
        root_layout.addLayout(location_layout)

        workspace_layout = QHBoxLayout()
        workspace_layout.addWidget(QLabel("Run directory"))
        self.workspace_edit = QLineEdit(str(UserDataPaths.discover().executions))
        workspace_layout.addWidget(self.workspace_edit, 1)
        workspace_browse = QPushButton("Browse...")
        workspace_browse.clicked.connect(self._choose_workspace)
        workspace_layout.addWidget(workspace_browse)
        root_layout.addLayout(workspace_layout)
        self.workspace_browse_button = workspace_browse
        self._execution_location_changed()

        self.case_table = QTableWidget(0, 4)
        self.case_table.setHorizontalHeaderLabels(
            ["Run", "Case", "Configuration", "Status"]
        )
        self.case_table.verticalHeader().setVisible(False)
        self.case_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.case_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.case_table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch
        )
        self.case_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.ResizeToContents
        )
        self.case_table.setColumnWidth(0, 58)
        self.case_table.setColumnWidth(3, 120)
        root_layout.addWidget(self.case_table, 3)

        action_layout = QHBoxLayout()
        select_all = QPushButton("Select all")
        select_all.clicked.connect(lambda: self._set_all_checks(True))
        select_none = QPushButton("Select none")
        select_none.clicked.connect(lambda: self._set_all_checks(False))
        refresh = QPushButton("Refresh cases")
        refresh.clicked.connect(self._load_cases)
        workflows = QPushButton("Workflows...")
        workflows.clicked.connect(self._open_workflows)
        mesh_layout_button = QPushButton("Mesh Layout...")
        mesh_layout_button.clicked.connect(self._open_mesh_layout)
        tests_button = QPushButton("Unit / Regression Tests...")
        tests_button.clicked.connect(self._open_tests)
        self.build_button = QPushButton("Build")
        self.build_button.clicked.connect(lambda: self._start_build(rebuild=False))
        self.rebuild_button = QPushButton("Rebuild")
        self.rebuild_button.clicked.connect(lambda: self._start_build(rebuild=True))
        history_button = QPushButton("Run History...")
        history_button.clicked.connect(self._open_run_history)
        self.run_button = QPushButton("Run")
        self.run_button.setToolTip("Run the selected catalog cases in isolated workspace copies.")
        self.run_button.setObjectName("primaryButton")
        self.run_button.clicked.connect(self._start_run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop_run)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        action_layout.addWidget(select_all)
        action_layout.addWidget(select_none)
        action_layout.addWidget(refresh)
        action_layout.addWidget(workflows)
        action_layout.addWidget(mesh_layout_button)
        action_layout.addWidget(self.build_button)
        action_layout.addWidget(self.rebuild_button)
        action_layout.addWidget(tests_button)
        action_layout.addWidget(history_button)
        action_layout.addStretch(1)
        action_layout.addWidget(self.run_button)
        action_layout.addWidget(self.stop_button)
        action_layout.addWidget(close_button)
        root_layout.addLayout(action_layout)

        root_layout.addWidget(QLabel("Execution log"))
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(10000)
        self.log_view.setPlaceholderText("Process output will appear here.")
        root_layout.addWidget(self.log_view, 2)

        self.setCentralWidget(central)
        self.statusBar().showMessage("Ready")
        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #f5f7fa; color: #202938; }
            QLabel#pageTitle { font-size: 24px; font-weight: 650; color: #172b4d; }
            QLabel#subtitle { color: #5d6b82; padding-bottom: 4px; }
            QLineEdit, QComboBox, QPlainTextEdit, QTableWidget {
                background: #ffffff; border: 1px solid #d5dce6; border-radius: 6px;
                padding: 6px;
            }
            QTableWidget { gridline-color: #e7ebf1; }
            QHeaderView::section {
                background: #edf1f6; color: #35445b; border: none;
                border-bottom: 1px solid #d5dce6; padding: 8px; font-weight: 600;
            }
            QPushButton {
                background: #ffffff; border: 1px solid #cbd4e1;
                border-radius: 6px; padding: 8px 12px;
            }
            QPushButton:hover { background: #eaf0f8; }
            QPushButton:disabled { color: #9aa4b2; background: #edf0f4; }
            QPushButton#primaryButton {
                background: #2459a9; color: white; border: 1px solid #2459a9;
                font-weight: 600; padding-left: 16px; padding-right: 16px;
            }
            QPushButton#primaryButton:hover { background: #1d4c93; }
            QStatusBar { background: #edf1f6; color: #526177; }
            """
        )

    def _refresh_executable_candidates(self) -> None:
        """Populate Debug and Release build slots plus existing installations."""
        current_text = self.executable_combo.currentText().strip()
        candidates = executable_candidates(self._repository_root)
        self.executable_combo.clear()
        self.executable_combo.addItems([str(path) for path in candidates])
        if current_text:
            self.executable_combo.setEditText(current_text)
        elif candidates:
            self.executable_combo.setCurrentIndex(0)

    def _load_cases(self) -> None:
        """Refresh the table from the shared catalog validator."""
        try:
            self._cases = load_cases(self._repository_root)
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.critical(self, "Cannot load cases", str(error))
            self.statusBar().showMessage("Case catalog validation failed")
            return

        self.case_table.blockSignals(True)
        self.case_table.setRowCount(0)
        self._row_by_id.clear()
        for case in self._cases:
            row = self.case_table.rowCount()
            self.case_table.insertRow(row)
            checkbox = QTableWidgetItem()
            checkbox.setFlags(checkbox.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            is_default_case = case.profile is not None
            is_selected = self._single_case_mode or is_default_case
            checkbox.setCheckState(
                Qt.CheckState.Checked if is_selected else Qt.CheckState.Unchecked
            )
            self.case_table.setItem(row, 0, checkbox)
            self.case_table.setItem(row, 1, QTableWidgetItem(case.title))
            config_text = "Typed profile" if case.profile else "Legacy config"
            self.case_table.setItem(row, 2, QTableWidgetItem(config_text))
            self.case_table.setItem(row, 3, QTableWidgetItem(RunState.PENDING.value))
            self._row_by_id[case.case_id] = row
        self.case_table.blockSignals(False)
        self.statusBar().showMessage(
            f"Loaded {len(self._cases)} case(s) from case_catalog.json"
        )

    def _set_all_checks(self, checked: bool) -> None:
        for row in range(self.case_table.rowCount()):
            item = self.case_table.item(row, 0)
            if item is not None:
                item.setCheckState(
                    Qt.CheckState.Checked
                    if checked
                    else Qt.CheckState.Unchecked
                )

    def _choose_executable(self) -> None:
        """Open the file picker in the currently entered path's directory."""
        current_text = self.executable_combo.currentText().strip()
        current_path = Path(current_text).expanduser() if current_text else None
        start_directory = (
            str(current_path.parent)
            if current_path is not None and current_path.parent.is_dir()
            else str(self._repository_root)
        )
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select OneFLOW executable",
            start_directory,
            "OneFLOW executable (OneFLOW.exe);;Executable (*.exe);;All files (*)",
        )
        if path:
            self.executable_combo.setEditText(path)

    def _execution_location_changed(self, *_args) -> None:
        custom = self.execution_location_combo.currentData() == "custom_workspace"
        self.workspace_edit.setEnabled(custom)
        self.workspace_browse_button.setEnabled(custom)

    def _choose_workspace(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select run directory", self.workspace_edit.text()
        )
        if path:
            self.workspace_edit.setText(path)

    def _start_run(self) -> None:
        if self._task_worker is not None and self._task_worker.isRunning():
            return
        if self._worker is not None and self._worker.isRunning():
            return
        executable_text = self.executable_combo.currentText().strip()
        if not executable_text or not Path(executable_text).expanduser().is_file():
            QMessageBox.warning(
                self,
                "Executable required",
                "Select a valid OneFLOW.exe from the path list or use Browse...",
            )
            return

        selected: list[CaseDefinition] = []
        for case in self._cases:
            row = self._row_by_id[case.case_id]
            item = self.case_table.item(row, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                selected.append(case)
        if not selected:
            QMessageBox.information(
                self, "No cases selected", "Select at least one case to run."
            )
            return

        execution_mode = self.execution_location_combo.currentData()
        if execution_mode == "custom_workspace":
            workspace_text = self.workspace_edit.text().strip()
            if not workspace_text:
                QMessageBox.warning(
                    self, "Run directory required", "Choose a writable run directory."
                )
                return
            workspace_root = Path(workspace_text).expanduser().resolve()
        elif execution_mode == "default_workspace":
            workspace_root = UserDataPaths.discover().executions
        else:
            QMessageBox.warning(
                self,
                "Invalid execution location",
                "Cases from the case catalog must run from a workspace copy.",
            )
            return

        try:
            workspace_root.relative_to(self._repository_root)
        except ValueError:
            pass
        else:
            QMessageBox.warning(
                self,
                "Workspace must be outside the source repository",
                "Run copies and generated files must stay outside the OneFLOW source "
                "repository. Choose a directory elsewhere.",
            )
            return

        try:
            workspace_root.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            QMessageBox.critical(self, "Workspace unavailable", str(error))
            return

        run_id = create_run_id()
        started_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        log_path = (UserDataPaths.discover().logs / f"{run_id}.log").resolve()
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(
                f"Case Workbench run {run_id}\nExecutable: {Path(executable_text).resolve()}\n",
                encoding="utf-8",
            )
            initial_record = WorkflowRunRecord(
                run_id=run_id,
                workflow_id="case-workbench",
                status="running",
                started_at=started_at,
                executable=str(Path(executable_text).resolve()),
                steps=({
                    "type": "case_queue",
                    "case_ids": [case.case_id for case in selected],
                    "execution_mode": execution_mode,
                    "workspace": str(workspace_root / run_id),
                },),
                message="Case queue started.",
                log_path=str(log_path),
            )
            self._history_repository.save_run(initial_record)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot record run", str(error))
            return
        self._active_case_run_record = initial_record
        self._active_case_run_steps = []
        self._active_case_run_log_path = log_path
        self.log_view.appendPlainText(
            f"\nRun {run_id}\nExecutable: {Path(executable_text).resolve()}\n"
            + f"Workspace: {workspace_root / run_id}"
        )
        self.run_button.setEnabled(False)
        self.build_button.setEnabled(False)
        self.rebuild_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        for case in selected:
            self._set_case_status(case.case_id, RunState.PENDING.value)

        self._worker = RunQueueWorker(
            selected, Path(executable_text), workspace_root, run_id,
            execution_mode=execution_mode,
        )
        self._worker.log_line.connect(self._append_case_run_log)
        self._worker.case_state.connect(self._set_case_status)
        self._worker.case_result.connect(self._show_result)
        self._worker.finished.connect(self._queue_finished)
        self._worker.start()
        self.statusBar().showMessage(f"Running {len(selected)} case(s)...")

    def _start_build(self, *, rebuild: bool) -> None:
        """Build the selected Debug or Release configuration from the case window."""
        if self._worker is not None and self._worker.isRunning():
            return
        if self._task_worker is not None and self._task_worker.isRunning():
            return
        executable = self.executable_combo.currentText().strip()
        build_paths = {
            configuration: str(
                (self._repository_root / "build" / "bin" / configuration / "OneFLOW.exe").resolve()
            )
            for configuration in ("Debug", "Release")
        }
        build_config = next(
            (
                configuration
                for configuration, path in build_paths.items()
                if executable.casefold() == path.casefold()
            ),
            "Debug",
        )
        command = build_cmake_command(
            self._repository_root, build_config, rebuild=rebuild
        )
        self.log_view.clear()
        self.run_button.setEnabled(False)
        self.build_button.setEnabled(False)
        self.rebuild_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._task_worker = TestCommandWorker([command])
        self._task_worker.output_line.connect(self.log_view.appendPlainText)
        self._task_worker.completed.connect(self._build_completed)
        self._task_worker.finished.connect(self._build_worker_finished)
        self._task_worker.start()
        self.statusBar().showMessage(f"{'Rebuilding' if rebuild else 'Building'} {build_config}...")

    def _build_completed(self, succeeded: bool, message: str) -> None:
        self.log_view.appendPlainText(message)
        self.statusBar().showMessage(message)

    def _build_worker_finished(self) -> None:
        self._task_worker = None
        self.run_button.setEnabled(True)
        self.build_button.setEnabled(True)
        self.rebuild_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def _open_mesh_layout(self) -> None:
        dialog = MeshLayoutDialog(self._repository_root, self)
        dialog.exec()

    def _open_run_history(self) -> None:
        dialog = RunHistoryDialog(UserDataPaths.discover().workflow_data, self)
        dialog.exec()

    def _open_tests(self) -> None:
        executable = self.executable_combo.currentText().strip()
        build_paths = {
            configuration: str(
                (self._repository_root / "build" / "bin" / configuration / "OneFLOW.exe").resolve()
            )
            for configuration in ("Debug", "Release")
        }
        build_config = next(
            (
                configuration
                for configuration, path in build_paths.items()
                if executable.casefold() == path.casefold()
            ),
            "Debug",
        )
        dialog = TestRunnerDialog(
            self._repository_root,
            self,
            build_config=build_config,
            executable=executable,
        )
        dialog.exec()

    def _open_workflows(self) -> None:
        dialog = WorkflowDialog(
            self._repository_root,
            self.executable_combo.currentText().strip(),
            self,
        )
        dialog.exec()

    def _stop_run(self) -> None:
        stopping = False
        if self._worker is not None:
            self._worker.request_stop()
            stopping = True
        if self._task_worker is not None:
            self._task_worker.request_stop()
            stopping = True
        if stopping:
            self.stop_button.setEnabled(False)
            self.statusBar().showMessage(
                "Stop requested; waiting for the active process to exit."
            )

    def _set_case_status(self, case_id: str, status: str) -> None:
        row = self._row_by_id.get(case_id)
        if row is not None:
            self.case_table.item(row, 3).setText(status)

    def _append_case_run_log(self, line: str) -> None:
        self.log_view.appendPlainText(line)
        if self._active_case_run_log_path is None:
            return
        try:
            with self._active_case_run_log_path.open("a", encoding="utf-8") as stream:
                stream.write(line + "\n")
        except OSError as error:
            self.statusBar().showMessage(f"Could not write run log: {error}")

    def _show_result(self, result: CaseRunResult) -> None:
        self.log_view.appendPlainText(
            f"[{result.case_id}] {result.message}\nWorkspace: {result.workspace}"
        )
        self._active_case_run_steps.append({
            "type": "case",
            "case_id": result.case_id,
            "status": result.state.value.lower(),
            "return_code": result.return_code,
            "workspace": str(result.workspace),
            "message": result.message,
        })

    def _queue_finished(self) -> None:
        record = self._active_case_run_record
        if record is not None:
            statuses = {step["status"] for step in self._active_case_run_steps}
            if "failed" in statuses:
                status = "failed"
            elif "cancelled" in statuses:
                status = "cancelled"
            else:
                status = "succeeded"
            updated_record = WorkflowRunRecord(
                run_id=record.run_id,
                workflow_id=record.workflow_id,
                status=status,
                started_at=record.started_at,
                finished_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                executable=record.executable,
                steps=tuple(self._active_case_run_steps),
                artifacts=record.artifacts,
                message=f"Case queue finished with status '{status}'.",
                note=record.note,
                log_path=record.log_path,
            )
            try:
                self._history_repository.save_run(updated_record)
            except (OSError, ValueError) as error:
                self.log_view.appendPlainText(f"Could not save run history: {error}")
            self._active_case_run_record = None
        self.run_button.setEnabled(True)
        self.build_button.setEnabled(True)
        self.rebuild_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.statusBar().showMessage("Run queue finished")
        self._worker = None

    def closeEvent(self, event) -> None:
        """Avoid destroying a running QThread while its child process is active."""
        queue_running = self._worker is not None and self._worker.isRunning()
        task_running = self._task_worker is not None and self._task_worker.isRunning()
        if queue_running or task_running:
            self._stop_run()
            self.statusBar().showMessage("Stopping the active task before closing.")
            event.ignore()
            return
        event.accept()


def run_desktop_app() -> int:
    """Create the application and show the activity selector first."""
    from .home_window import WorkbenchHomeWindow

    application = QApplication.instance() or QApplication([])
    application.setApplicationName("OneFLOW Workbench")
    window = WorkbenchHomeWindow()
    window.show()
    return application.exec()
