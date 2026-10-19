"""Home screen for choosing a OneFLOW Workbench activity."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .executable_discovery import executable_candidates
from .mesh_layout_dialog import MeshLayoutDialog
from .run_history_dialog import RunHistoryDialog
from .test_dialog import TestRunnerDialog
from .user_data_paths import UserDataPaths
from .workflow_dialog import WorkflowDialog


class WorkbenchHomeWindow(QMainWindow):
    """Route users to an activity instead of forcing Case Workbench on startup."""

    def __init__(self, repository_root: Path | None = None) -> None:
        super().__init__()
        self._repository_root = (
            repository_root or Path(__file__).resolve().parents[2]
        ).resolve()
        self._case_window = None
        self.setWindowTitle("OneFLOW Workbench")
        self.resize(900, 650)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        root = QVBoxLayout(central)
        root.setContentsMargins(28, 24, 28, 24)
        root.setSpacing(16)

        title = QLabel("OneFLOW Workbench")
        title.setObjectName("pageTitle")
        subtitle = QLabel(
            "Choose what you want to run. Cases, workflow definitions, and run history "
            "are separate parts of the workbench."
        )
        subtitle.setObjectName("subtitle")
        subtitle.setWordWrap(True)
        root.addWidget(title)
        root.addWidget(subtitle)

        executable_row = QHBoxLayout()
        executable_row.addWidget(QLabel("OneFLOW executable"))
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(42)
        self.executable_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        executable_row.addWidget(self.executable_combo, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._choose_executable)
        executable_row.addWidget(browse)
        root.addLayout(executable_row)

        grid = QGridLayout()
        grid.setSpacing(12)
        self._add_activity(
            grid, 0, 0, "Run an existing case",
            "Select a case from the shared catalog and run it through the Case Workbench.",
            self._open_case_workbench,
        )
        self._add_activity(
            grid, 0, 1, "Configure multi-case run",
            "Configure an ordered multi-directory run. Cases do not declare dependencies.",
            self._open_multi_case_workflow,
        )
        self._add_activity(
            grid, 1, 0, "Run a dependency workflow",
            "Define producer/consumer relationships, declared outputs, and execution order.",
            self._open_dependency_workflow,
        )
        self._add_activity(
            grid, 1, 1, "Run tests",
            "Run unit tests and optionally selected CFD or mesh regression suites.",
            self._open_tests,
        )
        self._add_activity(
            grid, 2, 0, "Run history",
            "Review previous runs, logs, status, and recorded artifacts.",
            self._open_history,
        )
        self._add_activity(
            grid, 2, 1, "Mesh Layout Editor",
            "Create or edit mesh layouts and their case templates.",
            self._open_mesh_layout,
        )
        root.addLayout(grid, 1)

        footer = QLabel(
            "Execution modes are intentionally distinct: one case, multiple independent "
            "cases, or a dependency-aware workflow."
        )
        footer.setObjectName("footer")
        footer.setWordWrap(True)
        root.addWidget(footer)
        self.setCentralWidget(central)
        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #f5f7fa; color: #202938; }
            QLabel#pageTitle { font-size: 26px; font-weight: 650; color: #172b4d; }
            QLabel#subtitle, QLabel#footer { color: #5d6b82; }
            QComboBox { background: white; border: 1px solid #d5dce6;
                border-radius: 6px; padding: 7px; }
            QPushButton {
                background: white; border: 1px solid #cbd4e1; border-radius: 8px;
                padding: 12px; text-align: left; font-weight: 600; min-height: 28px;
            }
            QPushButton:hover { background: #eaf0f8; border-color: #8fa9cf; }
            QLabel.activityDescription { color: #5d6b82; padding: 0 4px 4px 4px; }
            """
        )

    def _add_activity(
        self,
        grid: QGridLayout,
        row: int,
        column: int,
        title: str,
        description: str,
        callback,
    ) -> None:
        card = QWidget(self)
        layout = QVBoxLayout(card)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)
        button = QPushButton(title)
        button.setMinimumHeight(48)
        button.clicked.connect(callback)
        detail = QLabel(description)
        detail.setObjectName("activityDescription")
        detail.setProperty("class", "activityDescription")
        detail.setWordWrap(True)
        layout.addWidget(button)
        layout.addWidget(detail)
        grid.addWidget(card, row, column)

    def _choose_executable(self) -> None:
        current = self.executable_combo.currentText().strip()
        current_path = Path(current).expanduser() if current else None
        start = (
            str(current_path.parent)
            if current_path is not None and current_path.parent.is_dir()
            else str(self._repository_root)
        )
        path, _ = QFileDialog.getOpenFileName(
            self, "Select OneFLOW executable", start,
            "OneFLOW executable (OneFLOW.exe);;Executable (*.exe);;All files (*)",
        )
        if path:
            self.executable_combo.setEditText(path)

    def _executable(self) -> str:
        return self.executable_combo.currentText().strip()

    def _open_case_workbench(self) -> None:
        from .main_window import MainWindow

        self._case_window = MainWindow(self._repository_root)
        executable = self._executable()
        if executable:
            self._case_window.executable_combo.setEditText(executable)
        self._case_window.show()
        self._case_window.raise_()
        self._case_window.activateWindow()

    def _open_multi_case_workflow(self) -> None:
        dialog = WorkflowDialog(
            self._repository_root,
            self._executable(),
            self,
            initial_mode="multi_directory",
        )
        dialog.exec()

    def _open_dependency_workflow(self) -> None:
        dialog = WorkflowDialog(
            self._repository_root,
            self._executable(),
            self,
            initial_mode="case_graph",
        )
        dialog.exec()

    def _open_tests(self) -> None:
        TestRunnerDialog(self._repository_root, self).exec()

    def _open_history(self) -> None:
        RunHistoryDialog(UserDataPaths.discover().workflow_data, self).exec()

    def _open_mesh_layout(self) -> None:
        MeshLayoutDialog(self._repository_root, self).exec()
