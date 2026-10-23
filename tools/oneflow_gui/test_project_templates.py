"""Tests for project template discovery and creation."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.project_templates import (
    available_project_templates,
    create_project,
)


class ProjectTemplatesTest(unittest.TestCase):
    def test_lists_blank_and_existing_mesh_templates(self) -> None:
        template_ids = {item.template_id for item in available_project_templates()}
        self.assertIn("blank", template_ids)
        self.assertIn("mesh_workflow", template_ids)
        self.assertIn("minimal_rectangle", template_ids)

    def test_creates_blank_project_with_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "repository"
            (repository / "grid").mkdir(parents=True)
            target = root / "workspace"

            project = create_project(repository, "blank", "MyProject", target)

            self.assertEqual(project.root, (target / "MyProject").resolve())
            metadata = json.loads(
                (project.root / ".oneflow-project.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["format"], "oneflow.project")
            self.assertEqual(metadata["name"], "MyProject")
            self.assertTrue((project.root / "README.md").is_file())

    def test_copies_template_without_modifying_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "repository"
            source = repository / "grid" / "mesh_workflow"
            source.mkdir(parents=True)
            (source / "case.json").write_text('{"task": "Grid"}\n', encoding="utf-8")
            target = root / "workspace"

            project = create_project(repository, "mesh_workflow", "ConvertMesh", target)

            self.assertEqual(
                (project.root / "case.json").read_text(encoding="utf-8"),
                '{"task": "Grid"}\n',
            )
            self.assertEqual(
                (source / "case.json").read_text(encoding="utf-8"),
                '{"task": "Grid"}\n',
            )

    def test_rejects_existing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "repository"
            repository.mkdir()
            target = root / "workspace"
            (target / "Existing").mkdir(parents=True)

            with self.assertRaises(FileExistsError):
                create_project(repository, "blank", "Existing", target)

    def test_rejects_destination_inside_repository(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "repository"
            repository.mkdir()

            with self.assertRaises(ValueError):
                create_project(repository, "blank", "Unsafe", repository / "grid")

    def test_rejects_invalid_name_and_unknown_template(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "repository"
            repository.mkdir()
            target = root / "workspace"

            with self.assertRaises(ValueError):
                create_project(repository, "blank", "../Unsafe", target)
            with self.assertRaises(ValueError):
                create_project(repository, "unknown", "Project", target)


if __name__ == "__main__":
    unittest.main()
