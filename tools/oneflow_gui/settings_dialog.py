"""Dialog for configuring the local OneFLOW source checkout."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from .application_settings import ApplicationSettings


class SettingsDialog(QDialog):
    """Edit user-level application settings without changing run data paths."""

    def __init__(self, settings: ApplicationSettings, parent=None) -> None:
        super().__init__(parent)
        self._settings = settings
        self._repository_root: Path | None = None
        self.setWindowTitle("OneFLOW Basic Configuration")
        self.resize(760, 150)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("OneFLOW source root"))
        row = QHBoxLayout()
        self.repository_edit = QLineEdit(str(settings.repository_root))
        self.repository_edit.setPlaceholderText("Select the OneFLOW source directory")
        row.addWidget(self.repository_edit, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse)
        row.addWidget(browse)
        layout.addLayout(row)

        note = QLabel(
            "Choose an existing checkout directory containing CMakeLists.txt. "
            "The setting is stored in your user data directory."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        actions = QHBoxLayout()
        actions.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        save = QPushButton("Save")
        save.clicked.connect(self._save)
        actions.addWidget(cancel)
        actions.addWidget(save)
        layout.addLayout(actions)

    @property
    def repository_root(self) -> Path | None:
        return self._repository_root

    def _browse(self) -> None:
        current = self.repository_edit.text().strip()
        start = current if Path(current).is_dir() else str(Path.home())
        selected = QFileDialog.getExistingDirectory(
            self, "Select OneFLOW source root", start
        )
        if selected:
            self.repository_edit.setText(selected)

    def _save(self) -> None:
        value = self.repository_edit.text().strip()
        if not value:
            QMessageBox.warning(self, "Invalid source root", "Select a source directory.")
            return
        root = Path(value).expanduser()
        if not root.is_dir() or not (root / "CMakeLists.txt").is_file():
            QMessageBox.warning(
                self,
                "Invalid source root",
                "The selected directory must exist and contain CMakeLists.txt.",
            )
            return
        self._repository_root = root.resolve()
        self.accept()
