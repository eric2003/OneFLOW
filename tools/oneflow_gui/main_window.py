"""Qt desktop interface for selecting and running OneFLOW cases."""

from __future__ import annotations

import threading
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
    ) -> None:
        super().__init__()
        self._cases = cases
        self._executable = executable
        self._workspace_root = workspace_root
        self._run_id = run_id
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
            )
            self.case_result.emit(result)
            self.case_state.emit(case.case_id, result.state.value)


class MainWindow(QMainWindow):
    """Desktop shell; case loading and process execution live in services."""

    def __init__(self, repository_root: Path | None = None) -> None:
        super().__init__()
        self._repository_root = (
            repository_root or Path(__file__).resolve().parents[2]
        ).resolve()
        self._cases: list[CaseDefinition] = []
        self._worker: RunQueueWorker | None = None
        self._row_by_id: dict[str, int] = {}

        self.setWindowTitle("OneFLOW Case Runner")
        self.resize(1100, 760)
        self._build_ui()
        self._load_cases()

    def _build_ui(self) -> None:
        """Build a compact layout with configuration, queue, and log areas."""
        central = QWidget(self)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(20, 18, 20, 18)
        root_layout.setSpacing(12)

        heading = QLabel("OneFLOW Case Runner")
        heading.setObjectName("pageTitle")
        subtitle = QLabel(
            "Select cases, run them in isolated workspaces, and review process logs."
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

        workspace_layout = QHBoxLayout()
        workspace_layout.addWidget(QLabel("Run workspace"))
        self.workspace_edit = QLineEdit(str(Path.home() / "OneFLOWRuns"))
        workspace_layout.addWidget(self.workspace_edit, 1)
        workspace_browse = QPushButton("Browse...")
        workspace_browse.clicked.connect(self._choose_workspace)
        workspace_layout.addWidget(workspace_browse)
        root_layout.addLayout(workspace_layout)

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
        self.run_button = QPushButton("Run selected cases")
        self.run_button.setObjectName("primaryButton")
        self.run_button.clicked.connect(self._start_run)
        self.stop_button = QPushButton("Stop queue")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop_run)
        action_layout.addWidget(select_all)
        action_layout.addWidget(select_none)
        action_layout.addWidget(refresh)
        action_layout.addStretch(1)
        action_layout.addWidget(self.run_button)
        action_layout.addWidget(self.stop_button)
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
        """Populate the path dropdown with existing conventional installations."""
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

        self.case_table.setRowCount(0)
        self._row_by_id.clear()
        for case in self._cases:
            row = self.case_table.rowCount()
            self.case_table.insertRow(row)
            checkbox = QTableWidgetItem()
            checkbox.setFlags(checkbox.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            checkbox.setCheckState(
                Qt.CheckState.Checked
                if case.profile is not None
                else Qt.CheckState.Unchecked
            )
            self.case_table.setItem(row, 0, checkbox)
            self.case_table.setItem(row, 1, QTableWidgetItem(case.title))
            config_text = "Typed profile" if case.profile else "Legacy config"
            self.case_table.setItem(row, 2, QTableWidgetItem(config_text))
            self.case_table.setItem(row, 3, QTableWidgetItem(RunState.PENDING.value))
            self._row_by_id[case.case_id] = row
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

    def _choose_workspace(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select run workspace", self.workspace_edit.text()
        )
        if path:
            self.workspace_edit.setText(path)

    def _start_run(self) -> None:
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

        workspace_text = self.workspace_edit.text().strip()
        if not workspace_text:
            QMessageBox.warning(
                self, "Workspace required", "Choose a writable run workspace."
            )
            return

        workspace_root = Path(workspace_text).expanduser().resolve()
        try:
            workspace_root.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            QMessageBox.critical(self, "Workspace unavailable", str(error))
            return

        run_id = create_run_id()
        self.log_view.appendPlainText(
            f"\nRun {run_id}\nExecutable: {Path(executable_text).resolve()}\n"
            f"Workspace: {workspace_root / run_id}"
        )
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        for case in selected:
            self._set_case_status(case.case_id, RunState.PENDING.value)

        self._worker = RunQueueWorker(
            selected, Path(executable_text), workspace_root, run_id
        )
        self._worker.log_line.connect(self.log_view.appendPlainText)
        self._worker.case_state.connect(self._set_case_status)
        self._worker.case_result.connect(self._show_result)
        self._worker.finished.connect(self._queue_finished)
        self._worker.start()
        self.statusBar().showMessage(f"Running {len(selected)} case(s)...")

    def _stop_run(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)
            self.statusBar().showMessage(
                "Stop requested; waiting for active process to exit."
            )

    def _set_case_status(self, case_id: str, status: str) -> None:
        row = self._row_by_id.get(case_id)
        if row is not None:
            self.case_table.item(row, 3).setText(status)

    def _show_result(self, result: CaseRunResult) -> None:
        self.log_view.appendPlainText(
            f"[{result.case_id}] {result.message}\nWorkspace: {result.workspace}"
        )

    def _queue_finished(self) -> None:
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.statusBar().showMessage("Run queue finished")
        self._worker = None

    def closeEvent(self, event) -> None:
        """Avoid destroying a running QThread while its child process is active."""
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop()
            self.statusBar().showMessage("Stopping the active case before closing.")
            event.ignore()
            return
        event.accept()


def run_desktop_app() -> int:
    """Create the application and enter Qt's event loop."""
    application = QApplication.instance() or QApplication([])
    application.setApplicationName("OneFLOW Case Runner")
    window = MainWindow()
    window.show()
    return application.exec()
