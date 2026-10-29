"""Home screen for choosing a OneFLOW Workbench activity."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .application_settings import (
    ApplicationSettings,
    resolve_repository_root,
    workflow_editor_mode_for_run_mode,
)
from .configuration_dialog import ConfigurationDialog
from .mesh_layout_dialog import MeshLayoutDialog
from .geometry_contract_dialog import GeometryContractDialog
from .repository_test_inventory_dialog import RepositoryTestInventoryDialog
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
from .workspace_repository import load_workspace
from .workspace_workflow_adapter import build_workflow_definition_from_workspace
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
        self._test_runner_dialog = None
        self._project_windows = []
        self.setWindowTitle("OneFLOW Workbench")
        self.resize(900, 650)
        self._build_ui()
        QTimer.singleShot(0, self._review_active_background_tasks)

    def _build_ui(self) -> None:
        central = QWidget(self)
        root = QVBoxLayout(central)
        root.setContentsMargins(36, 32, 36, 32)
        root.setSpacing(20)

        title = QLabel("OneFLOW Workbench")
        title.setObjectName("pageTitle")
        root.addWidget(title)

        project_button = QPushButton("Project")
        project_button.setMinimumHeight(58)
        project_button.setMinimumWidth(240)
        project_button.setMaximumWidth(360)
        project_button.clicked.connect(self._open_project_menu)
        root.addWidget(project_button)

        for label, callback in (
            ("Tests...", self._open_tests),
            ("Repository Test Inventory...", self._open_test_inventory),
            ("Geometry Contracts...", self._open_geometry_contracts),
            ("Mesh Layout Editor...", self._open_mesh_layout),
            ("Background Tasks...", self._open_background_tasks),
            ("Configuration...", self._open_configuration),
        ):
            button = QPushButton(label)
            button.setMinimumHeight(44)
            button.setMaximumWidth(360)
            button.clicked.connect(callback)
            root.addWidget(button)

        root.addStretch(1)

        self.setCentralWidget(central)
        self.setStyleSheet(
            """
            QMainWindow, QWidget { background: #f5f7fa; color: #202938; }
            QLabel#pageTitle { font-size: 26px; font-weight: 650; color: #172b4d; }
            QPushButton {
                background: white; border: 1px solid #cbd4e1; border-radius: 8px;
                padding: 12px; text-align: left; font-weight: 600; min-height: 28px;
            }
            QPushButton:hover { background: #eaf0f8; border-color: #8fa9cf; }
            """
        )

    def _open_project_menu(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Project")
        dialog.setMinimumWidth(420)
        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(8)

        primary_actions = (
            ("New Project...", self._new_project),
            ("Open Project...", self._open_project),
            ("Open Multi-Case Workspace...", self._open_workspace),
            ("Create Multi-Case Workspace...", self._create_workspace),
            ("Run History...", self._open_history),
            ("Saved Workflows...", self._open_saved_workflows),
        )
        for label, callback in primary_actions:
            button = QPushButton(label)
            button.setMinimumHeight(42)
            button.clicked.connect(
                lambda checked=False, action=callback: self._run_project_action(
                    dialog, action
                )
            )
            layout.addWidget(button)

        dialog.setStyleSheet(
            """
            QDialog, QWidget { background: #f5f7fa; color: #202938; }
            QPushButton {
                background: white; border: 1px solid #cbd4e1; border-radius: 7px;
                padding: 10px; text-align: left; font-weight: 600;
            }
            QPushButton:hover { background: #eaf0f8; border-color: #8fa9cf; }
            """
        )
        dialog.exec()

    def _run_project_action(self, dialog: QDialog, callback) -> None:
        dialog.accept()
        QTimer.singleShot(0, callback)

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
        paths = UserDataPaths.discover()
        start_directory = paths.projects if paths.projects.is_dir() else Path.home()
        directory = QFileDialog.getExistingDirectory(
            self,
            "Open OneFLOW project directory",
            str(start_directory),
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
        paths = UserDataPaths.discover()
        start_directory = paths.workspaces if paths.workspaces.is_dir() else Path.home()
        path, _ = QFileDialog.getOpenFileName(
            self, "Open OneFLOW workspace", str(start_directory),
            "OneFLOW workspace (*.ofw)",
        )
        if path:
            self._edit_workspace(Path(path))

    def _create_workspace(self) -> None:
        self._edit_workspace()

    def _edit_workspace(self, workspace_path: Path | None = None) -> None:
        dialog = WorkspaceDialog(self, workspace_path)
        if dialog.exec() != QDialog.DialogCode.Accepted or not dialog.run_requested:
            return
        if dialog.workspace_path is None:
            return
        try:
            workspace = load_workspace(dialog.workspace_path)
            workflow = build_workflow_definition_from_workspace(
                workspace, dialog.workspace_path
            )
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot prepare workspace run", str(error))
            return
        WorkflowDialog(
            self._repository_root, self._executable(), self, initial_workflow=workflow
        ).exec()

    def _open_saved_workflows(self) -> None:
        WorkflowDialog(
            self._repository_root,
            self._executable(),
            self,
            initial_mode=workflow_editor_mode_for_run_mode(self._settings.run_mode),
        ).exec()

    def _open_tests(self) -> None:
        if self._test_runner_dialog is not None:
            if self._test_runner_dialog.isVisible():
                self._test_runner_dialog.raise_()
                self._test_runner_dialog.activateWindow()
                return
            self._test_runner_dialog = None

        dialog = TestRunnerDialog(
            self._repository_root,
            self,
            build_config=self._settings.build_configuration,
            executable=self._executable(),
        )
        self._test_runner_dialog = dialog
        dialog.finished.connect(self._test_runner_closed)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _test_runner_closed(self, _result: int) -> None:
        self._test_runner_dialog = None

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

    def _open_test_inventory(self) -> None:
        RepositoryTestInventoryDialog(self._repository_root, self).exec()

    def _open_geometry_contracts(self) -> None:
        GeometryContractDialog(self._repository_root, self).exec()

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
