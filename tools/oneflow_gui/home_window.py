"""Home screen for choosing a OneFLOW Workbench activity."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
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
from .executable_discovery import executable_candidates
from .mesh_layout_dialog import MeshLayoutDialog
from .project_repository import open_project
from .project_window import ProjectWindow
from .new_project_dialog import NewProjectDialog
from .run_history_dialog import RunHistoryDialog
from .test_dialog import TestRunnerDialog
from .settings_dialog import SettingsDialog
from .user_data_paths import UserDataPaths
from .workflow_dialog import WorkflowDialog
from .workspace_dialog import WorkspaceDialog


class WorkbenchHomeWindow(QMainWindow):
    """Route users to an activity instead of forcing Case Workbench on startup."""

    def __init__(self, repository_root: Path | None = None) -> None:
        super().__init__()
        discovered_root = repository_root or Path(__file__).resolve().parents[2]
        self._settings = ApplicationSettings.load()
        self._repository_root = resolve_repository_root(
            self._settings, discovered_root
        )
        self._case_window = None
        self._project_windows = []
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
            "Open or create a project workspace, choose how its cases execute, or run "
            "Workbench tasks such as catalog cases, unit tests, and CFD/grid regression."
        )
        subtitle.setObjectName("subtitle")
        subtitle.setWordWrap(True)
        root.addWidget(title)
        root.addWidget(subtitle)

        project_actions = QHBoxLayout()
        new_project_button = QPushButton("New Project...")
        new_project_button.setMinimumHeight(38)
        new_project_button.clicked.connect(self._new_project)
        project_actions.addWidget(new_project_button)
        project_actions.addStretch(1)
        root.addLayout(project_actions)

        configuration_row = QHBoxLayout()
        configuration_row.addWidget(QLabel("Build configuration"))
        self.build_config_combo = QComboBox()
        self.build_config_combo.addItems(["Debug", "Release"])
        self.build_config_combo.setToolTip(
            "Choose the build configuration used by build and unit-test tasks."
        )
        configuration_row.addWidget(self.build_config_combo, 1)
        root.addLayout(configuration_row)

        executable_row = QHBoxLayout()
        executable_row.addWidget(QLabel("OneFLOW executable"))
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setToolTip(
            "Choose the Debug or Release build for development, or browse to another executable."
        )
        self.executable_combo.setMinimumContentsLength(42)
        self.executable_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        self._build_executable_paths = {
            configuration: str(
                (self._repository_root / "build" / "bin" / configuration / "OneFLOW.exe").resolve()
            )
            for configuration in ("Debug", "Release")
        }
        self.build_config_combo.currentTextChanged.connect(
            self._build_configuration_changed
        )
        executable_row.addWidget(self.executable_combo, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._choose_executable)
        executable_row.addWidget(browse)
        root.addLayout(executable_row)

        grid = QGridLayout()
        grid.setSpacing(12)
        self._add_activity(
            grid, 0, 0, "Open Project",
            "Open any existing project directory, whether or not it belongs to a built-in catalog.",
            self._open_project,
        )
        self._add_activity(
            grid, 0, 1, "Case Database",
            "Run all selected catalog cases from protected workspace copies; all cases are selected by default.",
            self._open_case_database,
        )

        run_mode_row = QHBoxLayout()
        run_mode_row.addWidget(QLabel("Run mode"))
        self.run_mode_combo = QComboBox()
        self.run_mode_combo.addItem("Run single case", "single_case")
        self.run_mode_combo.addItem("Multi-case run", "multi_case")
        self.run_mode_combo.addItem("Dependency workflow", "dependency_workflow")
        self.run_mode_combo.setMinimumHeight(38)
        run_mode_row.addWidget(self.run_mode_combo, 1)
        open_run_mode = QPushButton("Continue...")
        open_run_mode.setMinimumHeight(38)
        open_run_mode.clicked.connect(self._open_selected_run_mode)
        run_mode_row.addWidget(open_run_mode)
        root.addLayout(run_mode_row)

        self._add_activity(
            grid, 1, 0, "Task Runner",
            "Run unit tests or CFD/grid regression; catalog cases are another task type.",
            self._open_tests,
        )
        self._add_activity(
            grid, 1, 1, "Run history",
            "Review previous runs, logs, status, and recorded artifacts.",
            self._open_history,
        )
        self._add_activity(
            grid, 2, 0, "Mesh Layout Editor",
            "Create or edit mesh layouts and their case templates.",
            self._open_mesh_layout,
        )
        self._add_activity(
            grid, 2, 1, "Saved workflows",
            "Open existing multi-case or dependency workflow definitions and their run history.",
            self._open_saved_workflows,
        )
        self._add_activity(
            grid, 3, 0, "Open Workspace",
            "Open a saved .ofw definition containing one or more cases and an execution policy.",
            self._open_workspace,
        )
        self._add_activity(
            grid, 3, 1, "Create Workspace",
            "Create a workspace and register case directories for single, sequential, independent, or dependency-based execution.",
            self._create_workspace,
        )
        self._add_activity(
            grid, 4, 0, "Basic Configuration",
            "Configure the OneFLOW source checkout used by build, test, and case tools.",
            self._open_settings,
        )
        root.addLayout(grid, 1)

        footer = QLabel(
            "A workspace can represent one case or a collection of cases. Execution policy belongs "
            "to the workspace. Catalog cases, unit tests, and CFD/grid regression are task types "
            "with different execution locations; their controls are being unified incrementally."
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

    def _build_configuration_changed(self, configuration: str) -> None:
        current = self.executable_combo.currentText().strip()
        if current.casefold() in {
            path.casefold() for path in self._build_executable_paths.values()
        }:
            self.executable_combo.setEditText(self._build_executable_paths[configuration])

    def _open_settings(self) -> None:
        dialog = SettingsDialog(self._settings, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        repository_root = dialog.repository_root
        if repository_root is None:
            return
        settings = ApplicationSettings(repository_root)
        try:
            settings.save()
        except OSError as error:
            QMessageBox.critical(self, "Cannot save settings", str(error))
            return
        self._settings = settings
        self._repository_root = repository_root
        self.executable_combo.clear()
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        self._build_configuration_changed(self.build_config_combo.currentText())
        QMessageBox.information(
            self,
            "Basic Configuration",
            f"OneFLOW source root set to:\n{self._repository_root}",
        )

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

    def _open_case_database(self) -> None:
        from .main_window import MainWindow

        self._case_window = MainWindow(
            self._repository_root,
            single_case_mode=True,
            build_config=self.build_config_combo.currentText(),
        )
        executable = self._executable()
        if executable:
            self._case_window.executable_combo.setEditText(executable)
        self._case_window.show()
        self._case_window.raise_()
        self._case_window.activateWindow()

    def _open_single_case(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select project for single-case run",
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
        executable = self._executable()
        if executable and hasattr(window, "executable_combo"):
            window.executable_combo.setEditText(executable)
        window.show()
        window.raise_()
        window.activateWindow()

    def _open_selected_run_mode(self) -> None:
        mode = self.run_mode_combo.currentData()
        if mode == "single_case":
            self._open_single_case()
        elif mode == "multi_case":
            WorkflowDialog(
                self._repository_root,
                self._executable(),
                self,
                initial_mode="multi_directory",
            ).exec()
        elif mode == "dependency_workflow":
            WorkflowDialog(
                self._repository_root,
                self._executable(),
                self,
                initial_mode="case_graph",
            ).exec()

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
            build_config=self.build_config_combo.currentText(),
            executable=self._executable(),
        ).exec()

    def _open_history(self) -> None:
        RunHistoryDialog(UserDataPaths.discover().workflow_data, self).exec()

    def _open_mesh_layout(self) -> None:
        MeshLayoutDialog(self._repository_root, self).exec()
