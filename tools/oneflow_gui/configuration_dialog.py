"""Unified configuration dialog for OneFLOW Workbench."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from .application_settings import ApplicationSettings
from .executable_discovery import executable_candidates


class ConfigurationDialog(QDialog):
    """Edit and persist the Workbench's shared run and build preferences."""

    def __init__(
        self,
        settings: ApplicationSettings,
        default_repository_root: Path,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("OneFLOW Configuration")
        self.setMinimumWidth(680)
        self._settings = settings
        self._default_repository_root = default_repository_root.resolve()
        self._saved_values = None

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Configure the source checkout, executable, run mode, and build configuration."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.repository_root_edit = QLineEdit(str(settings.repository_root))
        self.repository_root_edit.setToolTip(
            "Source checkout containing CMakeLists.txt."
        )
        root_row = QHBoxLayout()
        root_row.addWidget(self.repository_root_edit, 1)
        browse_root = QPushButton("Browse...")
        browse_root.clicked.connect(self._choose_repository_root)
        root_row.addWidget(browse_root)
        form.addRow("OneFLOW source root", root_row)

        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(42)
        self.executable_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(settings.repository_root)]
        )
        self.executable_combo.setEditText(
            settings.executable_path or self._default_executable(
                settings.repository_root, settings.build_configuration
            )
        )
        self._uses_automatic_executable = self._is_automatic_executable(
            settings.repository_root, self.executable_combo.currentText()
        )
        self.executable_combo.currentTextChanged.connect(
            self._executable_selection_changed
        )
        executable_row = QHBoxLayout()
        executable_row.addWidget(self.executable_combo, 1)
        browse_executable = QPushButton("Browse...")
        browse_executable.clicked.connect(self._choose_executable)
        executable_row.addWidget(browse_executable)
        form.addRow("OneFLOW executable", executable_row)

        self.run_mode_combo = QComboBox()
        self.run_mode_combo.addItem("Run single case", "single_case")
        self.run_mode_combo.addItem("Multi-case run", "multi_case")
        self.run_mode_combo.addItem("Dependency workflow", "dependency_workflow")
        self.run_mode_combo.setCurrentIndex(
            max(0, self.run_mode_combo.findData(settings.run_mode))
        )
        form.addRow("Run mode", self.run_mode_combo)

        self.build_config_combo = QComboBox()
        self.build_config_combo.addItems(["Debug", "Release"])
        self.build_config_combo.setCurrentText(settings.build_configuration)
        form.addRow("Build configuration", self.build_config_combo)
        layout.addLayout(form)

        self.build_config_combo.currentTextChanged.connect(
            self._build_configuration_changed
        )
        self.repository_root_edit.textChanged.connect(self._source_root_changed)
        self.repository_root_edit.editingFinished.connect(self._refresh_candidates_from_edit)

        actions = QHBoxLayout()
        self.reset_button = QPushButton("Reset to defaults")
        self.reset_button.clicked.connect(self._reset)
        actions.addWidget(self.reset_button)
        actions.addStretch(1)
        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.accepted.connect(self._save)
        self.buttons.rejected.connect(self.reject)
        actions.addWidget(self.buttons)
        layout.addLayout(actions)

        self._saved_values = self._current_values()

    @staticmethod
    def _default_executable(repository_root: Path, configuration: str) -> str:
        return str(
            (repository_root / "build" / "bin" / configuration / "OneFLOW.exe").resolve()
        )

    @property
    def settings(self) -> ApplicationSettings:
        return self._settings

    def _current_values(self) -> tuple[str, str, str, str]:
        return (
            self.repository_root_edit.text().strip(),
            self.executable_combo.currentText().strip(),
            self.run_mode_combo.currentData(),
            self.build_config_combo.currentText(),
        )

    def _choose_repository_root(self) -> None:
        current = self.repository_root_edit.text().strip()
        start = current if current and Path(current).is_dir() else str(
            self._default_repository_root
        )
        directory = QFileDialog.getExistingDirectory(
            self, "Select OneFLOW source root", start
        )
        if directory:
            self.repository_root_edit.setText(directory)

    def _source_root_changed(self, value: str) -> None:
        if self._uses_automatic_executable and value.strip():
            self.executable_combo.setEditText(
                self._default_executable(
                    Path(value.strip()), self.build_config_combo.currentText()
                )
            )

    def _is_automatic_executable(self, repository_root: Path, value: str) -> bool:
        return value.strip().casefold() in {
            self._default_executable(repository_root, configuration).casefold()
            for configuration in ("Debug", "Release")
        }

    def _executable_selection_changed(self, value: str) -> None:
        root_text = self.repository_root_edit.text().strip()
        if root_text:
            self._uses_automatic_executable = self._is_automatic_executable(
                Path(root_text), value
            )

    def _refresh_candidates_from_edit(self) -> None:
        value = self.repository_root_edit.text().strip()
        if value and Path(value).is_dir():
            self._refresh_executable_candidates(Path(value))

    def _refresh_executable_candidates(self, repository_root: Path) -> None:
        current = self.executable_combo.currentText().strip()
        self.executable_combo.clear()
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(repository_root)]
        )
        self.executable_combo.setEditText(current)

    def _choose_executable(self) -> None:
        current = self.executable_combo.currentText().strip()
        current_path = Path(current).expanduser() if current else None
        start = (
            str(current_path.parent)
            if current_path is not None and current_path.parent.is_dir()
            else self.repository_root_edit.text().strip()
        )
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select OneFLOW executable",
            start,
            "OneFLOW executable (OneFLOW.exe);;Executable (*.exe);;All files (*)",
        )
        if path:
            self.executable_combo.setEditText(path)

    def _build_configuration_changed(self, configuration: str) -> None:
        root_text = self.repository_root_edit.text().strip()
        if root_text and self._uses_automatic_executable:
            self.executable_combo.setEditText(
                self._default_executable(Path(root_text), configuration)
            )

    def _reset(self) -> None:
        root = self._default_repository_root
        self.repository_root_edit.setText(str(root))
        self.build_config_combo.setCurrentText("Debug")
        self.executable_combo.setEditText(self._default_executable(root, "Debug"))
        self.run_mode_combo.setCurrentIndex(
            self.run_mode_combo.findData("single_case")
        )
        self._refresh_executable_candidates(root)

    def _save(self, close_dialog: bool = True) -> bool:
        root_text = self.repository_root_edit.text().strip()
        root = Path(root_text).expanduser() if root_text else None
        if root is None or not root.is_dir() or not (root / "CMakeLists.txt").is_file():
            QMessageBox.warning(
                self,
                "Invalid source root",
                "Select an existing OneFLOW source directory containing CMakeLists.txt.",
            )
            return False
        root = root.resolve()
        executable = self.executable_combo.currentText().strip()
        if not executable:
            executable = self._default_executable(
                root, self.build_config_combo.currentText()
            )
        settings = ApplicationSettings(
            root,
            self.build_config_combo.currentText(),
            executable,
            self.run_mode_combo.currentData(),
        )
        try:
            settings.save()
        except OSError as error:
            QMessageBox.critical(self, "Cannot save configuration", str(error))
            return False
        self._settings = settings
        self._saved_values = self._current_values()
        if close_dialog:
            self.accept()
        return True

    def closeEvent(self, event) -> None:
        if self._current_values() == self._saved_values:
            event.accept()
            return
        prompt = QMessageBox(self)
        prompt.setWindowTitle("Unsaved configuration")
        prompt.setText("Configuration has unsaved changes.")
        save_button = prompt.addButton("Save", QMessageBox.ButtonRole.AcceptRole)
        discard_button = prompt.addButton(
            "Discard", QMessageBox.ButtonRole.DestructiveRole
        )
        cancel_button = prompt.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
        prompt.setDefaultButton(save_button)
        prompt.exec()
        clicked = prompt.clickedButton()
        if clicked is save_button:
            if self._save(close_dialog=False):
                event.accept()
            else:
                event.ignore()
        elif clicked is discard_button:
            event.accept()
        elif clicked is cancel_button:
            event.ignore()
        else:
            event.ignore()
