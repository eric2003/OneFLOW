"""Qt editor and runner for persisted multi-directory workflows."""

from __future__ import annotations

import json
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
from .workflow_graph_runner import run_case_graph
from .workflow_plan import build_case_graph_plan
from .workflow_process import run_process
from .user_data_paths import UserDataPaths
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
        *,
        case_graph: list[dict] | None = None,
        base_directory: Path | None = None,
        run_root: Path | None = None,
        log_path: Path | None = None,
    ) -> None:
        super().__init__()
        self._executable = executable
        self._directories = directories
        self._workflow_id = workflow_id
        self._run_id = run_id
        self._case_graph = case_graph
        self._base_directory = base_directory
        self._run_root = run_root
        self._log_path = log_path
        self._cancel = threading.Event()

    def request_stop(self) -> None:
        self._cancel.set()

    def run(self) -> None:
        started_at = _timestamp()
        command: list[str] = []
        return_code: int | None = None
        message = ""
        status = "failed"
        step_records: list[dict] = []
        artifacts: list[str] = []
        try:
            if not self._executable.is_file():
                raise FileNotFoundError(f"OneFLOW executable not found: {self._executable}")
            if self._case_graph is not None:
                if self._base_directory is None or self._run_root is None:
                    raise ValueError("Case graph execution requires a base directory and run directory.")
                plan = build_case_graph_plan(self._case_graph, self._base_directory)
                process_records: dict[str, dict] = {}

                def run_case(step):
                    case_command = [str(self._executable), "r", str(step.directory)]
                    self.log_line.emit(f"Running case '{step.case_id}'")
                    self.log_line.emit("Command arguments:")
                    for index, argument in enumerate(case_command):
                        self.log_line.emit(f"  [{index}] {argument}")
                    result = run_process(
                        case_command,
                        cwd=self._executable.parent,
                        on_output=self.log_line.emit,
                        is_cancelled=self._cancel.is_set,
                    )
                    process_records[step.case_id] = {
                        "command": case_command,
                        "return_code": result.return_code,
                    }
                    return result.status, result.message

                graph_result = run_case_graph(
                    plan,
                    run_case,
                    run_root=self._run_root,
                    is_cancelled=self._cancel.is_set,
                )
                status = graph_result.status
                for case_result in graph_result.cases:
                    record = {
                        "type": "case",
                        "case_id": case_result.case_id,
                        "status": case_result.status,
                        "message": case_result.message,
                    }
                    record.update(process_records.get(case_result.case_id, {}))
                    step_records.append(record)
                    self.log_line.emit(
                        f"Case '{case_result.case_id}': {case_result.status} - {case_result.message}"
                    )
                for plan_step in plan.steps:
                    staged_directory = self._run_root / "cases" / plan_step.case_id
                    for output in plan_step.outputs:
                        output_path = staged_directory / output.relative_path
                        if output_path.is_file():
                            artifacts.append(str(output_path))
                message = f"Case graph finished with status '{status}'."
            else:
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
                step_records = [
                    {
                        "type": "multi_directory",
                        "directories": [str(path) for path in self._directories],
                        "command": command,
                        "return_code": return_code,
                    }
                ]
        except (OSError, ValueError, KeyError, TypeError) as error:
            message = str(error)
        self.completed.emit(
            {
                "run_id": self._run_id,
                "workflow_id": self._workflow_id,
                "status": status,
                "started_at": started_at,
                "finished_at": _timestamp(),
                "executable": str(self._executable),
                "steps": step_records,
                "artifacts": artifacts,
                "message": message,
                "log_path": str(self._log_path) if self._log_path is not None else None,
            }
        )


class WorkflowDialog(QDialog):
    """Manage workflow definitions and persistent execution history."""

    def __init__(
        self,
        repository_root: Path,
        executable: str = "",
        parent: QWidget | None = None,
        initial_mode: str | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("OneFLOW Workflows")
        self.resize(900, 720)
        self._repository_root = repository_root
        self._user_data_paths = UserDataPaths.discover()
        self._repository = WorkflowRepository(
            self._user_data_paths.workflow_data
        )
        self._worker: WorkflowRunWorker | None = None
        self._active_log_path: Path | None = None
        self._selected_id: str | None = None
        self._build_ui(executable)
        self._reload_workflows()
        if initial_mode in {"multi_directory", "case_graph"}:
            self._new_workflow()
            self.mode_combo.setCurrentIndex(self.mode_combo.findData(initial_mode))

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
        self.mode_combo = QComboBox()
        self.mode_combo.addItem("Multi-directory", "multi_directory")
        self.mode_combo.addItem("Dependency graph", "case_graph")
        self.mode_combo.currentIndexChanged.connect(self._workflow_mode_changed)
        self.directory_edit = QPlainTextEdit()
        self.directory_edit.setPlaceholderText(
            "One directory per line. Order is preserved in the OneFLOW command."
        )
        self.directory_edit.setMinimumHeight(120)
        directory_row = QHBoxLayout()
        directory_row.addWidget(self.directory_edit, 1)
        self.directory_browse = QPushButton("Add directory...")
        self.directory_browse.clicked.connect(self._add_workflow_directory)
        directory_row.addWidget(self.directory_browse)
        directory_widget = QWidget()
        directory_widget.setLayout(directory_row)
        form.addRow("Workflow ID", self.id_edit)
        form.addRow("Name", self.name_edit)
        self.form_layout = form
        self.graph_edit = QPlainTextEdit()
        self.graph_edit.setPlaceholderText("JSON array of case objects, including dependencies, outputs, and inputs.")
        self.graph_edit.setMinimumHeight(160)
        self.base_directory_edit = QLineEdit()
        self.base_directory_edit.setPlaceholderText("Base directory containing the configured case directories")
        base_directory_row = QHBoxLayout()
        base_directory_row.addWidget(self.base_directory_edit, 1)
        self.base_directory_browse = QPushButton("Browse...")
        self.base_directory_browse.clicked.connect(self._choose_base_directory)
        base_directory_row.addWidget(self.base_directory_browse)
        base_directory_widget = QWidget()
        base_directory_widget.setLayout(base_directory_row)
        form.addRow("Description", self.description_edit)
        form.addRow("Workflow type", self.mode_combo)
        form.addRow("Directories", directory_widget)
        form.addRow("Case graph JSON", self.graph_edit)
        form.addRow("Graph base directory", base_directory_widget)
        root.addLayout(form)
        self._workflow_mode_changed()

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


    def _add_workflow_directory(self) -> None:
        """Append a directory selected from the filesystem to the workflow."""
        lines = self.directory_edit.toPlainText().splitlines()
        start = Path(lines[-1]).expanduser() if lines and lines[-1].strip() else self._repository_root
        if start.is_file():
            start = start.parent
        if not start.is_dir():
            start = self._repository_root
        directory = QFileDialog.getExistingDirectory(
            self, "Add workflow directory", str(start)
        )
        if directory:
            current = self.directory_edit.toPlainText().rstrip()
            self.directory_edit.setPlainText(f"{current}\n{directory}" if current else directory)


    def _choose_base_directory(self) -> None:
        """Choose the root directory used to resolve case graph paths."""
        current = self.base_directory_edit.text().strip()
        start = Path(current).expanduser() if current else self._repository_root
        if start.is_file():
            start = start.parent
        if not start.is_dir():
            start = self._repository_root
        directory = QFileDialog.getExistingDirectory(
            self, "Select workflow base directory", str(start)
        )
        if directory:
            self.base_directory_edit.setText(directory)

    def _workflow_mode_changed(self, *_args) -> None:
        graph_mode = self.mode_combo.currentData() == "case_graph"
        self.directory_edit.setVisible(not graph_mode)
        self.graph_edit.setVisible(graph_mode)
        self.base_directory_edit.setVisible(graph_mode)
        directory_label = self.form_layout.labelForField(self.directory_edit)
        graph_label = self.form_layout.labelForField(self.graph_edit)
        base_label = self.form_layout.labelForField(self.base_directory_edit)
        if directory_label is not None:
            directory_label.setVisible(not graph_mode)
        if graph_label is not None:
            graph_label.setVisible(graph_mode)
        if base_label is not None:
            base_label.setVisible(graph_mode)

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
            step_types = [step.get("type") for step in workflow.steps]
            if step_types == ["case_graph"]:
                graph_step = workflow.steps[0]
                self.mode_combo.setCurrentIndex(self.mode_combo.findData("case_graph"))
                self.graph_edit.setPlainText(json.dumps(graph_step.get("cases", []), indent=2))
                self.base_directory_edit.setText(str(graph_step.get("base_directory", "")))
                self.directory_edit.clear()
            else:
                directories = read_multi_directory_step(workflow.steps)
                self.mode_combo.setCurrentIndex(self.mode_combo.findData("multi_directory"))
                self.directory_edit.setPlainText("\n".join(directories))
                self.graph_edit.clear()
                self.base_directory_edit.clear()
        except (ValueError, TypeError) as error:
            self.directory_edit.clear()
            self.graph_edit.clear()
            self.log_view.appendPlainText(
                f"Cannot load workflow '{workflow.workflow_id}': {error}"
            )
        self.delete_button.setEnabled(True)
        self._reload_history(workflow.workflow_id)

    def _new_workflow(self) -> None:
        self.workflow_list.clearSelection()
        self.workflow_list.setCurrentRow(-1)
        self._selected_id = None
        self.id_edit.clear()
        self.id_edit.setReadOnly(False)
        self.name_edit.clear()
        self.description_edit.clear()
        self.mode_combo.setCurrentIndex(self.mode_combo.findData("multi_directory"))
        self.directory_edit.clear()
        self.graph_edit.clear()
        self.base_directory_edit.clear()
        self.delete_button.setEnabled(False)
        self.history_list.clear()
        self.name_edit.setFocus()

    def _save_workflow(self) -> bool:
        name = self.name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Name required", "Enter a workflow name.")
            return False
        workflow_id = self.id_edit.text().strip() or _new_id(name)
        if self.mode_combo.currentData() == "case_graph":
            try:
                graph_data = json.loads(self.graph_edit.toPlainText())
                imported_base = ""
                if isinstance(graph_data, list):
                    cases = graph_data
                elif isinstance(graph_data, dict) and isinstance(graph_data.get("steps"), list):
                    graph_steps = graph_data["steps"]
                    if len(graph_steps) != 1 or graph_steps[0].get("type") != "case_graph":
                        raise ValueError("Imported workflow JSON must contain exactly one case_graph step.")
                    cases = graph_steps[0].get("cases")
                    imported_base = str(graph_steps[0].get("base_directory", ""))
                elif isinstance(graph_data, dict) and graph_data.get("type") == "case_graph":
                    cases = graph_data.get("cases")
                    imported_base = str(graph_data.get("base_directory", ""))
                else:
                    raise ValueError("Case graph JSON must be a case array or a case_graph workflow definition.")
                if not isinstance(cases, list) or not cases:
                    raise ValueError("Case graph JSON must contain a non-empty cases array.")
                graph_step = {"type": "case_graph", "cases": cases}
                base_directory = self.base_directory_edit.text().strip() or imported_base
                if base_directory:
                    graph_step["base_directory"] = base_directory
                steps = (graph_step,)
            except (json.JSONDecodeError, ValueError, TypeError) as error:
                QMessageBox.warning(self, "Invalid case graph", str(error))
                return False
        else:
            directories = [line.strip() for line in self.directory_edit.toPlainText().splitlines() if line.strip()]
            if not directories:
                QMessageBox.warning(self, "Directories required", "Enter at least one directory.")
                return False
            steps = ({"type": "multi_directory", "directories": directories},)
        definition = WorkflowDefinition(
            workflow_id=workflow_id,
            name=name,
            description=self.description_edit.text().strip(),
            steps=steps,
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
        case_graph = None
        base_directory = None
        run_root = None
        directories: list[Path] = []
        try:
            workflow = self._repository.load_workflow(self._selected_id)
            executable_text = self.executable_edit.text().strip()
            if not executable_text:
                raise ValueError("Select a OneFLOW executable.")
            executable = Path(executable_text).expanduser().resolve()
            if not executable.is_file():
                raise FileNotFoundError(f"OneFLOW executable not found: {executable}")
            if workflow.steps[0].get("type") == "case_graph":
                graph_step = workflow.steps[0]
                case_graph = graph_step["cases"]
                base_text = str(graph_step.get("base_directory", "")).strip()
                if not base_text:
                    base_text = self.base_directory_edit.text().strip()
                if not base_text:
                    raise ValueError("Select a base directory for the case graph.")
                base_directory = Path(base_text).expanduser().resolve()
                if not base_directory.is_dir():
                    raise NotADirectoryError(f"Workflow base directory not found: {base_directory}")
                plan = build_case_graph_plan(case_graph, base_directory)
                missing = [step.directory for step in plan.steps if not step.directory.is_dir()]
                if missing:
                    raise NotADirectoryError(
                        "These workflow case directories do not exist:\n" +
                        "\n".join(str(path) for path in missing)
                    )
                command = [str(executable), "r", "<each staged case directory>"]
            else:
                directories = [
                    Path(value).expanduser().resolve()
                    for value in read_multi_directory_step(workflow.steps)
                ]
                missing = [path for path in directories if not path.is_dir()]
                if missing:
                    raise NotADirectoryError(
                        "These workflow directories do not exist:\n" +
                        "\n".join(str(path) for path in missing)
                    )
                command = build_multi_directory_command(executable, directories)
        except (OSError, ValueError, KeyError, StopIteration, TypeError) as error:
            QMessageBox.warning(self, "Invalid workflow", str(error))
            return
        run_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid4().hex[:8]
        started_at = _timestamp()
        log_path = (self._user_data_paths.logs / f"{run_id}.log").resolve()
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(f"Run {run_id}: {workflow.name}\n", encoding="utf-8")
        self._active_log_path = log_path
        if case_graph is not None:
            run_root = (self._user_data_paths.executions / run_id).resolve()
            initial_steps = ({
                "type": "case_graph",
                "base_directory": str(base_directory),
                "cases": case_graph,
                "command_template": command,
            },)
        else:
            initial_steps = ({
                "type": "multi_directory",
                "directories": [str(path) for path in directories],
                "command": command,
            },)
        self._repository.save_run(
            WorkflowRunRecord(
                run_id=run_id,
                workflow_id=workflow.workflow_id,
                status="running",
                started_at=started_at,
                executable=str(executable),
                steps=initial_steps,
                message="Workflow process started.",
                log_path=str(log_path),
            )
        )
        self.log_view.appendPlainText(f"\nRun {run_id}: {workflow.name}")
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = WorkflowRunWorker(
            executable,
            directories,
            workflow.workflow_id,
            run_id,
            case_graph=case_graph,
            base_directory=base_directory,
            run_root=run_root,
            log_path=log_path,
        )
        self._worker.log_line.connect(self._append_run_log)
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

    def _append_run_log(self, line: str) -> None:
        self.log_view.appendPlainText(line)
        if self._active_log_path is None:
            return
        try:
            with self._active_log_path.open("a", encoding="utf-8", newline="\n") as stream:
                stream.write(line + "\n")
        except OSError as error:
            self.log_view.appendPlainText(f"Could not write run log: {error}")

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
