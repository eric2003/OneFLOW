"""Tests for user-owned mesh project creation and reopening."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.mesh_project import (
    create_mesh_project_from_stl,
    create_mesh_project_from_template,
    import_mesh_project_directory,
    load_mesh_project,
    mesh_project_templates,
    replace_project_stl,
    save_mesh_project_settings,
)


class MeshProjectTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.home = Path(self.temporary.name)
        self.repository_root = Path(__file__).resolve().parents[2]
        self.environment = patch.dict(os.environ, {"ONEFLOW_HOME": str(self.home)})
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_prolate_spheroid_is_a_selectable_project_template(self) -> None:
        choices = mesh_project_templates(self.repository_root)
        self.assertIn(
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            {case_id for case_id, _title in choices},
        )
        self.assertIn(
            "stl_to_gmsh_python_sphere_automatic_tetrahedra",
            {case_id for case_id, _title in choices},
        )

    def test_create_and_reopen_project_from_prolate_template(self) -> None:
        project = create_mesh_project_from_template(
            self.repository_root,
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            "prolate_spheroid",
        )
        self.assertTrue(project["stl_path"].is_file())
        self.assertEqual(project["directory"].parent, self.home / "mesh-generation")
        self.assertEqual(project["mesh_case_id"], "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra")
        reopened = load_mesh_project(project["directory"])
        self.assertEqual(reopened["output_path"], project["output_path"])
        saved = save_mesh_project_settings(
            reopened,
            {"mesh_size": 0.12, "boundary_name": "body_wall"},
        )
        manifest = json.loads((saved["directory"] / "project.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["mesh_size"], 0.12)
        self.assertEqual(manifest["boundary_name"], "body_wall")

    def test_replacing_project_stl_updates_private_geometry_copy(self) -> None:
        source = self.repository_root / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python/input/model.stl"
        alternate = self.home / "alternate.stl"
        alternate.write_text("solid replacement\\nendsolid replacement\\n", encoding="ascii")
        project = create_mesh_project_from_template(
            self.repository_root,
            "stl_to_gmsh_python_sphere_automatic_tetrahedra",
            "replaceable",
        )
        updated = replace_project_stl(project, alternate)
        self.assertEqual(updated["stl_path"].read_text(encoding="ascii"), alternate.read_text(encoding="ascii"))
        self.assertNotEqual(updated["stl_path"], alternate)
        self.assertEqual(updated["mesh_case_id"], project["mesh_case_id"])

    def test_import_prepared_project_directory_and_reopen_it(self) -> None:
        source = self.home / "prepared_case"
        (source / "input").mkdir(parents=True)
        (source / "script").mkdir()
        stl = self.repository_root / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python/input/model.stl"
        (source / "input/model.stl").write_bytes(stl.read_bytes())
        (source / "script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3, "mesh_size": 0.15, "boundary_name": "body", "surface_name": "fluid"}),
            encoding="utf-8",
        )
        project = import_mesh_project_directory(self.repository_root, source)
        self.assertTrue(project["stl_path"].is_file())
        self.assertNotEqual(project["stl_path"], source / "input/model.stl")
        self.assertEqual(project["mesh_size"], 0.15)
        self.assertEqual(project["boundary_name"], "body")
        self.assertEqual(load_mesh_project(project["directory"])["output_path"], project["output_path"])

    def test_import_stl_copies_database_geometry_into_independent_project(self) -> None:
        source = self.repository_root / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python/input/model.stl"
        project = create_mesh_project_from_stl(self.repository_root, source, "imported_spheroid")
        self.assertTrue(project["stl_path"].is_file())
        self.assertNotEqual(project["stl_path"], source)
        self.assertEqual(project["mesh_case_id"], "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra")


if __name__ == "__main__":
    unittest.main()
