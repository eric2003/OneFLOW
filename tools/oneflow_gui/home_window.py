"""Home screen for choosing a OneFLOW Workbench activity."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QDialog,
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

from .application_settings import ApplicationSettings, resolve_repository_root
from .configuration_dialog import ConfigurationDialog
from .mesh_layout_dialog import MeshLayoutDialog
from .project_repository import open_project
from .project_window import ProjectWindow
from .new_project_dialog import NewProjectDialog
from .run_history_dialog import RunHistoryDialog
from .background_task_dialog import BackgroundTaskDialog
from .background_task_registry import BackgroundTaskRegistry
from .background_task_manager import request_background_stop
from .test_dialog import TestRunnerDialog
from .user_data_paths import UserDataPaths
from .workflow_dialog import WorkflowDialog
from .workspace_dialog import WorkspaceDialog


class WorkbenchHomeWindow(QMainWindow):
    """Route users to an activity instead of forcing Case Workbench on startup."""

    def __init__(self, repository_root: Path | None = None) -> None:
        super().__init__()
        discovered_root = repository_root or Path(__file__).resolve().parents[2]
        self._default_repository_root = discovered_root.resolve()
        self._settings = ApplicationSettings.load()
        self._repository_root = resolve_repository_root(
            self._settings, discovered_root
        )
        if not self._settings.executable_path:
            self._settings = ApplicationSettings(
                self._repository_root,
                self._settings.build_configuration,
                str(
                    (self._repository_root / "build" / "bin" / self._settings.build_configuration / "OneFLOW.exe").resolve()
                ),
                self._settings.run_mode,
            )
        self._case_window = None
        self._project_windows = []
        self.setWindowTitle("OneFLOW Workbench")
        self.resize(900, 650)
        self._build_ui()
        QTimer.singleShot(0, self._review_active_background_tasks)

    def _build_ui(self) -> None:
        central = QWidget(self)
        root = QVBoxLayout(central)
        root.setContentsMargins(28, 24, 28, 24)
        root.setSpacing(16)

        title = QLabel("OneFLOW Workbench")
        title.setObjectName("pageTitle")
        subtitle = QLabel(
            "Create or open a project folder to edit and run calculations. Multi-Case Workspace "
            "files (.ofw) remain available for existing case collections; other Workbench tasks "
            "include catalog cases, unit tests, and CFD/grid regression."
        )
        subtitle.setObjectName("subtitle")
        subtitle.setWordWrap(True)
        root.addWidget(title)
        root.addWidget(subtitle)

        project_actions = QHBoxLayout()
        new_project_button = QPushButton("New Project...")
        new_project_button.setMinimumHeight(48)
        new_project_button.setMinimumWidth(180)
        new_project_button.clicked.connect(self._new_project)
        project_actions.addWidget(new_project_button)

        open_project_button = QPushButton("Open Project...")
        open_project_button.setMinimumHeight(48)
        open_project_button.setMinimumWidth(180)
        open_project_button.clicked.connect(self._open_project)
        project_actions.addWidget(open_project_button)

        self.background_tasks_button = QPushButton("Background Tasks")
        self.background_tasks_button.setMinimumHeight(48)
        self.background_tasks_button.setMinimumWidth(180)
        self.background_tasks_button.setToolTip(
            "Manage active calculations. Finished task records can be shown from the manager."
        )
        self.background_tasks_button.clicked.connect(self._open_background_tasks)
        project_actions.addWidget(self.background_tasks_button)

        configuration_button = QPushButton("Configuration...")
        configuration_button.setMinimumHeight(48)
        configuration_button.setMinimumWidth(160)
        configuration_button.setToolTip(
            "Configure the source root, executable, run mode, and build configuration."
        )
        configuration_button.clicked.connect(self._open_configuration)
        project_actions.addWidget(configuration_button)

        project_actions.addStretch(1)
        root.addLayout(project_actions)

        grid = QGridLayout()
        grid.setSpacing(12)
        self._add_activity(
            grid, 0, 0, "Tests",
            "Run unit tests, CFD/grid regression, or selected cases from the Case Database.",
            self._open_tests,
        )
        self._add_activity(
            grid, 0, 1, "Run history",
            "Review previous runs, logs, status, and recorded artifacts.",
            self._open_history,
        )
        self._add_activity(
            grid, 1, 0, "Mesh Layout Editor",
            "Create or edit mesh layouts and their case templates.",
            self._open_mesh_layout,
        )
        self._add_activity(
            grid, 1, 1, "Saved workflows",
            "Open existing multi-case or dependency workflow definitions and their run history.",
            self._open_saved_workflows,
        )
        self._add_activity(
            grid, 2, 0, "Open Multi-Case Workspace...",
            "Open a saved .ofw case collection and its execution policy. Project folders use Open Project.",
            self._open_workspace,
        )
        self._add_activity(
            grid, 2, 1, "Create Multi-Case Workspace...",
            "Create a .ofw case collection for existing case directories. For a new templated calculation project, use New Project above.",
            self._create_workspace,
        )
        root.addLayout(grid, 1)

        footer = QLabel(
            "A Project is a calculation directory that can be created from a template or opened directly. A Multi-Case Workspace is a saved .ofw definition for grouping existing cases and selecting an execution policy. These are distinct legacy paths; new project workflows are the preferred path for new work."
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

    def _open_configuration(self) -> None:
        dialog = ConfigurationDialog(
            ApplicationSettings(
                self._repository_root,
                self._settings.build_configuration,
                self._settings.executable_path,
                self._settings.run_mode,
            ),
            self._default_repository_root,
            self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._settings = dialog.settings
        self._repository_root = self._settings.repository_root.resolve()
        self.statusBar().showMessage("Configuration saved.", 4000)

    def _executable(self) -> str:
        return self._settings.executable_path.strip()

    def _new_project(self) -> None:
        dialog = NewProjectDialog(self._repository_root, self)
        if dialog.exec() != QDialog.DialogCode.Accepted or dialog.project is None:
            return
        window = ProjectWindow(dialog.project, self._repository_root, self)
        self._project_windows.append(window)
        window.show()
        window.raise_()
        window.activateWindow()

    def _open_project(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self,
            "Open OneFLOW project directory",
            str(Path.home()),
        )
        if not directory:
            return
        try:
            project = open_project(directory)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot open project", str(error))
            return
        window = ProjectWindow(project, self._repository_root, self)
        self._project_windows.append(window)
        window.show()
        window.raise_()
        window.activateWindow()

    def _open_workspace(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open OneFLOW workspace", str(Path.home()),
            "OneFLOW workspace (*.ofw)",
        )
        if path:
            WorkspaceDialog(self, Path(path)).exec()

    def _create_workspace(self) -> None:
        WorkspaceDialog(self).exec()

    def _open_saved_workflows(self) -> None:
        WorkflowDialog(
            self._repository_root, self._executable(), self
        ).exec()

    def _open_tests(self) -> None:
        TestRunnerDialog(
            self._repository_root,
            self,
            build_config=self._settings.build_configuration,
            executable=self._executable(),
        ).exec()

    def _open_history(self) -> None:
        RunHistoryDialog(UserDataPaths.discover().workflow_data, self).exec()

    def _review_active_background_tasks(self) -> None:
        registry = BackgroundTaskRegistry(UserDataPaths.discover().root)
        active_tasks = registry.list_active_tasks()
        if not active_tasks:
            return
        self._open_background_tasks()

    def _open_background_tasks(self) -> None:
        BackgroundTaskDialog(UserDataPaths.discover().root, self).exec()

    def _open_mesh_layout(self) -> None:
        MeshLayoutDialog(self._repository_root, self).exec()

    def closeEvent(self, event) -> None:
        registry = BackgroundTaskRegistry(UserDataPaths.discover().root)
        active_tasks = registry.list_active_tasks()
        if not active_tasks:
            event.accept()
            return

        prompt = QMessageBox(self)
        prompt.setWindowTitle("OneFLOW calculations are running")
        prompt.setText(f"{len(active_tasks)} background task(s) are still running.")
        prompt.setInformativeText(
            "Choose whether to stop these calculations or leave them running after OneFLOW Workbench closes."
        )
        stop_button = prompt.addButton("Stop calculations and close", QMessageBox.ButtonRole.DestructiveRole)
        continue_button = prompt.addButton("Continue in background", QMessageBox.ButtonRole.AcceptRole)
        cancel_button = prompt.addButton("Keep Workbench open", QMessageBox.ButtonRole.RejectRole)
        prompt.setDefaultButton(stop_button)
        prompt.exec()
        clicked = prompt.clickedButton()
        if clicked is cancel_button:
            event.ignore()
            return
        if clicked is continue_button:
            for window in self._project_windows:
                window.set_application_shutdown_decision("continue")
            event.accept()
            return
        if clicked is not stop_button:
            event.ignore()
            return

        try:
            for task in active_tasks:
                request_background_stop(registry.root, task.task_id)
            for window in self._project_windows:
                window.set_application_shutdown_decision("stop")
        except (OSError, ValueError) as error:
            QMessageBox.critical(
                self, "Could not stop all calculations",
                f"Workbench will remain open so you can manage the tasks.\n{error}",
            )
            event.ignore()
            return
        event.accept()
