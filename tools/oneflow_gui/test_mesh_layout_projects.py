"""Tests for copying mesh layout template projects into user-owned workspaces."""

from __future__ import annotations

import json
import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.mesh_layout_projects import ensure_template_project


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


class MeshLayoutProjectsTest(unittest.TestCase):
    def test_template_copy_includes_project_files_and_preserves_source(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "projects"
            source_json = REPOSITORY_ROOT / "grid" / "turbulent_flat_plate" / "flat_plate.json"
            original = source_json.read_text(encoding="utf-8")

            project = ensure_template_project(REPOSITORY_ROOT, workspace, "flat_plate")

            self.assertEqual(project, (workspace / "turbulent_flat_plate").resolve())
            self.assertTrue((project / "flat_plate.json").is_file())
            self.assertTrue((project / "case.json").is_file())
            self.assertTrue((project / "script" / "control.txt").is_file())
            self.assertTrue((project / "grid" / "turbulentFlatPlateLayout.txt").is_file())
            self.assertEqual(source_json.read_text(encoding="utf-8"), original)

            edited = json.loads((project / "flat_plate.json").read_text(encoding="utf-8"))
            edited["plateLength"] = 3.0
            (project / "flat_plate.json").write_text(json.dumps(edited), encoding="utf-8")
            self.assertEqual(source_json.read_text(encoding="utf-8"), original)

    def test_existing_project_is_reused_without_overwriting_user_edits(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "projects"
            project = workspace / "minimal_rectangle"
            project.mkdir(parents=True)
            specification = project / "rectangle.json"
            specification.write_text('{"type":"rectangle","userEdit":true}', encoding="utf-8")

            result = ensure_template_project(REPOSITORY_ROOT, workspace, "rectangle")

            self.assertEqual(result, project.resolve())
            self.assertEqual(
                specification.read_text(encoding="utf-8"),
                '{"type":"rectangle","userEdit":true}',
            )
            self.assertFalse((project / "script").exists())

    def test_project_destination_inside_repository_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "outside the OneFLOW source tree"):
                ensure_template_project(
                    REPOSITORY_ROOT,
                    REPOSITORY_ROOT / "temporary-test-workspace",
                    "rectangle",
                )

    def test_laminar_flat_plate_template_is_available(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project = ensure_template_project(
                REPOSITORY_ROOT,
                Path(directory),
                "laminar_flat_plate",
            )
            self.assertTrue((project / "flat_plate.json").is_file())
            spec = json.loads((project / "flat_plate.json").read_text(encoding="utf-8"))
            self.assertEqual(spec["distributions"]["vertical"]["nearWallPointCount"], 41)

    def test_unknown_template_type_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError):
                ensure_template_project(REPOSITORY_ROOT, Path(directory), "unknown")


if __name__ == "__main__":
    unittest.main()
