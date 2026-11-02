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
from .user_data_paths import UserDataPaths


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
            "Create one project folder from a template under the selected parent folder. "
            "This does not create a separate .ofw multi-case workspace. Repository templates "
            "are copied; the source tree is never used as the new project's destination."
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
        self._template_default_name = ""
        self.template_combo.currentIndexChanged.connect(
            self._update_default_project_name
        )
        self._update_default_project_name()
        form.addRow("Project name", self.name_edit)

        destination_row = QHBoxLayout()
        self.destination_edit = QLineEdit(str(UserDataPaths.discover().projects))
        destination_row.addWidget(self.destination_edit, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._choose_destination)
        destination_row.addWidget(browse)
        form.addRow("Project parent folder", destination_row)
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

    def _update_default_project_name(self) -> None:
        template_id = self.template_combo.currentData()
        default_name = "LaminarFlatPlate" if template_id == "laminar_flat_plate" else ""
        current_name = self.name_edit.text().strip()
        if not current_name or current_name == self._template_default_name:
            self.name_edit.setText(default_name)
        self._template_default_name = default_name

    def _choose_destination(self) -> None
        path = QFileDialog.getExistingDirectory(
            self,
            "Choose parent directory for the new project",
            self.destination_edit.text() or str(Path.home()),
        )
        if path:
            self.destination_edit.setText(path)

    def _create_project(self) -> None:
        template_id = self.template_combo.currentData()
        project_name = self.name_edit.text()
        parent_directory = Path(self.destination_edit.text()).expanduser()
        destination = parent_directory / project_name.strip()
        overwrite = False

        if destination.exists() or destination.is_symlink():
            answer = QMessageBox.question(
                self,
                "Replace existing project?",
                "The project folder already exists:\n\n"
                + str(destination)
                + "\n\nReplacing it will permanently delete the entire existing "
                "project folder and then copy the template into a new folder. "
                "All existing files and results in that project will be lost.\n\n"
                "Do you want to replace it?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
            overwrite = True

        try:
            self.project = create_project(
                self._repository_root,
                template_id,
                project_name,
                parent_directory,
                overwrite=overwrite,
            )
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Cannot create project", str(error))
            return
        self.accept()
