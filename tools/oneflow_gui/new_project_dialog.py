"""Dialog for creating a new OneFLOW project from a template."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)
from PySide6.QtCore import Qt

from .project_templates import available_project_templates, create_project


class NewProjectDialog(QDialog):
    """Collect template and destination choices, then create a project."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self.project = None
        self.setWindowTitle("New OneFLOW Project")
        self.resize(620, 260)

        layout = QVBoxLayout(self)
        intro = QLabel(
            "Create a project in your workspace. Repository templates are copied; "
            "the source tree is never used as the new project's destination."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        form = QFormLayout()
        self.template_combo = QComboBox()
        for template in available_project_templates():
            self.template_combo.addItem(template.title, template.template_id)
            self.template_combo.setItemData(
                self.template_combo.count() - 1,
                template.description,
                Qt.ItemDataRole.ToolTipRole,
            )
        form.addRow("Project template", self.template_combo)

        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("e.g. my_mesh_project")
        form.addRow("Project name", self.name_edit)

        destination_row = QHBoxLayout()
        self.destination_edit = QLineEdit(str(Path.home() / "OneFLOWWorkspace"))
        destination_row.addWidget(self.destination_edit, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._choose_destination)
        destination_row.addWidget(browse)
        form.addRow("Parent directory", destination_row)
        layout.addLayout(form)

        actions = QHBoxLayout()
        actions.addStretch(1)
        create_button = QPushButton("Create Project")
        create_button.setDefault(True)
        create_button.clicked.connect(self._create_project)
        actions.addWidget(create_button)
        cancel_button = QPushButton("Cancel")
        cancel_button.clicked.connect(self.reject)
        actions.addWidget(cancel_button)
        layout.addLayout(actions)

    def _choose_destination(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self,
            "Choose parent directory for the new project",
            self.destination_edit.text() or str(Path.home()),
        )
        if path:
            self.destination_edit.setText(path)

    def _create_project(self) -> None:
        try:
            self.project = create_project(
                self._repository_root,
                self.template_combo.currentData(),
                self.name_edit.text(),
                self.destination_edit.text(),
            )
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Cannot create project", str(error))
            return
        self.accept()
