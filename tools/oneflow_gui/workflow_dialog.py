"""Qt editor and runner for persisted multi-directory workflows."""

from __future__ import annotations

import threading
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .workflow_execution import (
    build_multi_directory_command,
    read_multi_directory_step,
)
from .workflow_process import run_process
from .workflow_repository import (
    WorkflowDefinition,
    WorkflowRepository,
    WorkflowRunRecord,
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _new_id(name: str) -> str:
    slug = "".join(char.lower() if char.isascii() and char.isalnum() else "-" for char in name)
    slug = "-".join(part for part in slug.split("-") if part) or "workflow"
    return f"{slug[:40]}-{uuid4().hex[:8]}"


class WorkflowRunWorker(QThread):
    """Run one configured directory sequence in a single OneFLOW process."""

    log_line = Signal(str)
    completed = Signal(object)

    def __init__(
        self,
        executable: Path,
        directories: list[Path],
        workflow_id: str,
        run_id: str,
    ) -> None:
        super().__init__()
        self._executable = executable
        self._directories = directories
        self._workflow_id = workflow_id
        self._run_id = run_id
        self._cancel = threading.Event()

    def request_stop(self) -> None:
        self._cancel.set()

    def run(self) -> None:
        started_at = _timestamp()
        command: list[str] = []
        return_code: int | None = None
        message = ""
        status = "failed"
        try:
            if not self._executable.is_file():
                raise FileNotFoundError(f"OneFLOW executable not found: {self._executable}")
            for index, directory in enumerate(self._directories, start=1):
                if not directory.is_dir():
                    raise NotADirectoryError(f"Directory {index} does not exist: {directory}")
            command = build_multi_directory_command(self._executable, self._directories)
            self.log_line.emit("Command arguments:")
            for index, argument in enumerate(command):
                self.log_line.emit(f"  [{index}] {argument}")
            result = run_process(
                command,
                cwd=self._executable.parent,
                on_output=self.log_line.emit,
                is_cancelled=self._cancel.is_set,
            )
            return_code = result.return_code
            status = result.status
            message = result.message
        except (OSError, ValueError) as error:
            message = str(error)
        self.completed.emit(
            {
                "run_id": self._run_id,
                "workflow_id": self._workflow_id,
                "status": status,
                "started_at": started_at,
                "finished_at": _timestamp(),
                "executable": str(self._executable),
                "steps": [
                    {
                        "type": "multi_directory",
                        "directories": [str(path) for path in self._directories],
                        "command": command,
                        "return_code": return_code,
                    }
                ],
                "artifacts": [],
                "message": message,
            }
        )


class WorkflowDialog(QDialog):
    """Manage workflow definitions and persistent execution history."""

    def __init__(
        self,
        repository_root: Path,
        executable: str = "",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("OneFLOW Workflows")
        self.resize(900, 720)
        self._repository_root = repository_root
        self._repository = WorkflowRepository(
            Path.home() / ".oneflow" / "workflow-data"
        )
        self._worker: WorkflowRunWorker | None = None
        self._selected_id: str | None = None
        self._build_ui(executable)
        self._reload_workflows()

    def _build_ui(self, executable: str) -> None:
        root = QVBoxLayout(self)
        root.addWidget(QLabel("Saved workflows"))
        self.workflow_list = QListWidget()
        self.workflow_list.currentRowChanged.connect(self._select_workflow)
        root.addWidget(self.workflow_list, 1)

        form = QFormLayout()
        self.id_edit = QLineEdit()
        self.id_edit.setPlaceholderText("Generated for new workflows; stable after save")
        self.name_edit = QLineEdit()
        self.description_edit = QLineEdit()
        self.directory_edit = QPlainTextEdit()
        self.directory_edit.setPlaceholderText(
            "One directory per line. Order is preserved in the OneFLOW command."
        )
        self.directory_edit.setMinimumHeight(120)
        form.addRow("Workflow ID", self.id_edit)
        form.addRow("Name", self.name_edit)
        form.addRow("Description", self.description_edit)
        form.addRow("Directories", self.directory_edit)
        root.addLayout(form)

        executable_row = QHBoxLayout()
        executable_row.addWidget(QLabel("OneFLOW executable"))
        self.executable_edit = QLineEdit(executable)
        executable_row.addWidget(self.executable_edit, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._choose_executable)
        executable_row.addWidget(browse)
        root.addLayout(executable_row)

        buttons = QHBoxLayout()
        self.new_button = QPushButton("New")
        self.save_button = QPushButton("Save")
        self.delete_button = QPushButton("Delete")
        self.run_button = QPushButton("Run workflow")
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.new_button.clicked.connect(self._new_workflow)
        self.save_button.clicked.connect(self._save_workflow)
        self.delete_button.clicked.connect(self._delete_workflow)
        self.run_button.clicked.connect(self._start_workflow)
        self.stop_button.clicked.connect(self._stop_workflow)
        for button in (
            self.new_button, self.save_button, self.delete_button,
            self.run_button, self.stop_button,
        ):
            buttons.addWidget(button)
        root.addLayout(buttons)

        root.addWidget(QLabel("Persistent run history"))
        self.history_list = QListWidget()
        self.history_list.setMaximumHeight(140)
        root.addWidget(self.history_list)
        root.addWidget(QLabel("Execution log"))
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(10000)
        root.addWidget(self.log_view, 2)

    def _reload_workflows(self, select_id: str | None = None) -> None:
        workflows = self._repository.list_workflows()
        self.workflow_list.blockSignals(True)
        self.workflow_list.clear()
        for workflow in workflows:
            self.workflow_list.addItem(f"{workflow.name} [{workflow.workflow_id}]")
        self.workflow_list.blockSignals(False)
        if select_id:
            for index, workflow in enumerate(workflows):
                if workflow.workflow_id == select_id:
                    self.workflow_list.setCurrentRow(index)
                    return
        if workflows:
            self.workflow_list.setCurrentRow(0)
        else:
            self._new_workflow()

    def _select_workflow(self, row: int) -> None:
        workflows = self._repository.list_workflows()
        if row < 0 or row >= len(workflows):
            return
        workflow = workflows[row]
        self._selected_id = workflow.workflow_id
        self.id_edit.setText(workflow.workflow_id)
        self.id_edit.setReadOnly(True)
        self.name_edit.setText(workflow.name)
        self.description_edit.setText(workflow.description)
        try:
            directories = read_multi_directory_step(workflow.steps)
        except ValueError as error:
            directories = []
            self.log_view.appendPlainText(
                f"Cannot load workflow '{workflow.workflow_id}': {error}"
            )
        self.directory_edit.setPlainText("\\n".join(directories))
        self.delete_button.setEnabled(True)
        self._reload_history(workflow.workflow_id)

    def _new_workflow(self) -> None:
        self.workflow_list.clearSelection()
        self._selected_id = None
        self.id_edit.clear()
        self.id_edit.setReadOnly(False)
        self.name_edit.clear()
        self.description_edit.clear()
        self.directory_edit.clear()
        self.delete_button.setEnabled(False)
        self.history_list.clear()
        self.name_edit.setFocus()

    def _save_workflow(self) -> bool:
        name = self.name_edit.text().strip()
        directories = [line.strip() for line in self.directory_edit.toPlainText().splitlines() if line.strip()]
        if not name:
            QMessageBox.warning(self, "Name required", "Enter a workflow name.")
            return False
        if not directories:
            QMessageBox.warning(self, "Directories required", "Enter at least one directory.")
            return False
        workflow_id = self.id_edit.text().strip() or _new_id(name)
        definition = WorkflowDefinition(
            workflow_id=workflow_id,
            name=name,
            description=self.description_edit.text().strip(),
            steps=({"type": "multi_directory", "directories": directories},),
        )
        try:
            self._repository.save_workflow(definition)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot save workflow", str(error))
            return False
        self._selected_id = workflow_id
        self.id_edit.setText(workflow_id)
        self.id_edit.setReadOnly(True)
        self._reload_workflows(workflow_id)
        return True

    def _delete_workflow(self) -> None:
        workflow_id = self._selected_id
        if not workflow_id:
            return
        answer = QMessageBox.question(
            self,
            "Delete workflow",
            f"Delete '{workflow_id}'? Existing run history will be retained.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            self._repository.delete_workflow(workflow_id)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot delete workflow", str(error))
            return
        self._reload_workflows()

    def _choose_executable(self) -> None:
        current = Path(self.executable_edit.text()).expanduser()
        start = current.parent if current.parent.is_dir() else self._repository_root
        path, _ = QFileDialog.getOpenFileName(
            self, "Select OneFLOW executable", str(start),
            "OneFLOW executable (OneFLOW.exe);;Executable (*.exe);;All files (*)",
        )
        if path:
            self.executable_edit.setText(path)

    def _start_workflow(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        if not self._save_workflow():
            return
        if not self._selected_id:
            return
        try:
            workflow = self._repository.load_workflow(self._selected_id)
            directories = [
                Path(value).expanduser().resolve()
                for value in read_multi_directory_step(workflow.steps)
            ]
            executable_text = self.executable_edit.text().strip()
            if not executable_text:
                raise ValueError("Select a OneFLOW executable.")
            executable = Path(executable_text).expanduser().resolve()
            if not executable.is_file():
                raise FileNotFoundError(f"OneFLOW executable not found: {executable}")
            command = build_multi_directory_command(executable, directories)
        except (OSError, ValueError, KeyError, StopIteration, TypeError) as error:
            QMessageBox.warning(self, "Invalid workflow", str(error))
            return
        missing = [path for path in directories if not path.is_dir()]
        if missing:
            QMessageBox.warning(
                self, "Directory not found",
                "These workflow directories do not exist:\n" + "\n".join(str(path) for path in missing),
            )
            return
        run_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid4().hex[:8]
        started_at = _timestamp()
        self._repository.save_run(
            WorkflowRunRecord(
                run_id=run_id,
                workflow_id=workflow.workflow_id,
                status="running",
                started_at=started_at,
                executable=str(executable),
                steps=({"type": "multi_directory", "directories": [str(path) for path in directories], "command": command},),
                message="Workflow process started.",
            )
        )
        self.log_view.appendPlainText(f"\nRun {run_id}: {workflow.name}")
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = WorkflowRunWorker(executable, directories, workflow.workflow_id, run_id)
        self._worker.log_line.connect(self.log_view.appendPlainText)
        self._worker.completed.connect(self._workflow_completed)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()
        self._reload_history(workflow.workflow_id)

    def _stop_workflow(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)
            self.log_view.appendPlainText("Stop requested.")

    def _workflow_completed(self, data: dict) -> None:
        try:
            self._repository.save_run(WorkflowRunRecord.from_dict(data))
        except (OSError, ValueError) as error:
            self.log_view.appendPlainText(f"Could not persist run result: {error}")
        self.log_view.appendPlainText(data["message"])
        self._reload_history(data["workflow_id"])

    def _worker_finished(self) -> None:
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self._worker = None

    def _reload_history(self, workflow_id: str) -> None:
        self.history_list.clear()
        try:
            records = self._repository.list_runs(workflow_id)
        except (OSError, ValueError) as error:
            self.history_list.addItem(f"History unavailable: {error}")
            return
        for record in records:
            self.history_list.addItem(
                f"{record.started_at} | {record.status} | {record.run_id} | {record.message}"
            )

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop()
            QMessageBox.information(
                self, "Stopping workflow",
                "The active process is being stopped. Close the workflow window after it exits.",
            )
            event.ignore()
            return
        event.accept()
