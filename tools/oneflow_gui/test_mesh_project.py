"""Tests for user-owned mesh project creation and reopening."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.path_utils import paths_equal
from grid.stl_to_gmsh.pipe_bend_geometry import generate as generate_pipe_bend
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

    def test_internal_flow_cylinder_is_a_selectable_project_template(self) -> None:
        template_id = "stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra"
        choices = {case_id for case_id, _title in mesh_project_templates(self.repository_root)}
        self.assertIn(template_id, choices)
        project = create_mesh_project_from_template(
            self.repository_root, template_id, "cylinder_internal_flow"
        )
        self.assertTrue(project["stl_path"].is_file())
        self.assertEqual(project["mesh_case_id"], template_id)
        copied_manifest = (
            project["directory"] / "template/input/model.geometry.json"
        )
        self.assertTrue(copied_manifest.is_file())
        metadata = json.loads(copied_manifest.read_text(encoding="utf-8"))
        self.assertEqual(metadata["asset_id"], "geometry.cylinder.internal_flow")
        self.assertTrue(metadata["surface"]["watertight"])
        self.assertTrue(
            paths_equal(
                project["directory"].parent,
                self.home / "mesh-generation",
            )
        )

    def test_create_and_reopen_project_from_prolate_template(self) -> None:
        project = create_mesh_project_from_template(
            self.repository_root,
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            "prolate_spheroid",
        )
        self.assertTrue(project["stl_path"].is_file())
        self.assertTrue(paths_equal(project["directory"].parent, self.home / "mesh-generation"))
        self.assertEqual(project["mesh_case_id"], "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra")
        reopened = load_mesh_project(project["directory"])
        self.assertTrue(paths_equal(reopened["output_path"], project["output_path"]))
        saved = save_mesh_project_settings(
            reopened,
            {"mesh_size": 0.12, "boundary_name": "body_wall"},
        )
        manifest = json.loads((saved["directory"] / "project.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["mesh_size"], 0.12)
        self.assertEqual(manifest["boundary_name"], "body_wall")

    def test_replacing_project_stl_updates_private_geometry_copy(self) -> None:
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
        self.assertIsNone(updated["template_case_id"])
        self.assertEqual(
            updated["mesh_case_id"], "stl_to_gmsh_generic_3d_mesh_handoff"
        )

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
        self.assertTrue(paths_equal(load_mesh_project(project["directory"])["output_path"], project["output_path"]))


    def test_cylinder_region_manifest_is_preserved_in_mesh_project(self) -> None:
        cylinder = self.home / "cylinder_internal_flow.stl"
        cylinder.write_bytes(b"placeholder")
        sidecar = cylinder.with_suffix(".geometry.json")
        sidecar.write_text(
            json.dumps({
                "schema": "oneflow.stl-geometry-manifest",
                "version": 1,
                "geometry": "cylinder_internal_flow",
                "patch_triangle_indices": {
                    "wall": [0, 1], "inflow": [2], "outflow": [3]
                },
            }),
            encoding="utf-8",
        )
        project = create_mesh_project_from_stl(
            self.repository_root, cylinder, "internal_flow"
        )
        self.assertEqual(project["geometry"], "cylinder_internal_flow")
        copied_sidecar = project["stl_path"].with_suffix(".geometry.json")
        self.assertTrue(copied_sidecar.is_file())
        copied_manifest = json.loads(copied_sidecar.read_text(encoding="utf-8"))
        self.assertEqual(
            copied_manifest["patch_triangle_indices"],
            {"wall": [0, 1], "inflow": [2], "outflow": [3]},
        )

        updated = replace_project_stl(project, cylinder)
        self.assertEqual(updated["geometry"], "cylinder_internal_flow")
        self.assertTrue(updated["stl_path"].with_suffix(".geometry.json").is_file())


    def test_pipe_bend_region_manifest_and_patch_files_are_preserved(self) -> None:
        source_directory = self.home / "generated_pipe_bend"
        source_directory.mkdir()
        source_stl = source_directory / "pipe_bend_90_degree.stl"
        generated = generate_pipe_bend(source_stl)
        project = create_mesh_project_from_stl(
            self.repository_root, source_stl, "pipe_bend_90_degree"
        )
        self.assertEqual(project["geometry"], "pipe_bend_90_degree")
        self.assertEqual(project["mesh_case_id"], "stl_to_gmsh_pipe_bend_90_degree")
        self.assertEqual(project["boundary_regions"], ["inflow", "outflow", "wall"])
        copied_stl = project["stl_path"]
        copied_manifest = copied_stl.parent / "region_manifest.json"
        self.assertTrue(copied_manifest.is_file())
        region_manifest = json.loads(copied_manifest.read_text(encoding="utf-8"))
        self.assertEqual(
            region_manifest["closed_surface_sha256"],
            hashlib.sha256(copied_stl.read_bytes()).hexdigest(),
        )
        for name, record in region_manifest["regions"].items():
            region_file = copied_stl.parent / record["file"]
            self.assertTrue(region_file.is_file(), name)
            self.assertEqual(
                hashlib.sha256(region_file.read_bytes()).hexdigest(),
                record["sha256"],
            )
        copied_geometry_manifest = json.loads(
            copied_stl.with_suffix(".geometry.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            copied_geometry_manifest["region_manifest"], "region_manifest.json"
        )
        self.assertTrue(generated["region_manifest_path"].is_file())

    def test_project_reopen_recovers_verified_geometry_and_validation_lineage(self) -> None:
        source = self.home / "generated.stl"
        source.write_bytes(b"project-stl")
        source_digest = hashlib.sha256(source.read_bytes()).hexdigest()
        source.with_suffix(".geometry.json").write_text(
            json.dumps({
                "schema": "oneflow.stl-geometry-manifest",
                "asset_id": "geometry.sphere.unit",
                "stl_sha256": source_digest,
                "definition_sha256": "a" * 64,
                "generator": {"name": "oneflow_stl_geometry", "version": "1.0"},
            }),
            encoding="utf-8",
        )
        project = create_mesh_project_from_stl(
            self.repository_root, source, "lineage_project"
        )
        project["output_path"].write_bytes(b"project-msh")
        msh_digest = hashlib.sha256(project["output_path"].read_bytes()).hexdigest()
        validation_root = project["directory"] / "validation"
        msh_dir = validation_root / "stl_to_gmsh/results"
        msh_dir.mkdir(parents=True)
        staged_msh = msh_dir / project["output_path"].name
        staged_msh.write_bytes(project["output_path"].read_bytes())
        cgns_dir = validation_root / "gmsh_to_cgns/results"
        cgns_dir.mkdir(parents=True)
        cgns = cgns_dir / "lineage.cgns"
        cgns.write_bytes(b"cgns-artifact")
        cgns_digest = hashlib.sha256(cgns.read_bytes()).hexdigest()
        geometry_lineage = {
            **project["geometry_lineage"],
            "derived_mesh": {
                "format": "Gmsh MSH", "sha256": msh_digest,
                "parent": "source_representation",
            },
        }
        (cgns_dir / (cgns.name + ".geometry-lineage.json")).write_text(
            json.dumps({
                "schema": "oneflow.geometry-artifact-lineage",
                "version": 1,
                "artifact": {
                    "role": "CGNS volume mesh", "file_name": cgns.name,
                    "format": "CGNS", "sha256": cgns_digest,
                },
                "geometry_lineage": geometry_lineage,
                "parent_artifact": {
                    "role": "Gmsh MSH", "file_name": staged_msh.name,
                    "path": staged_msh.relative_to(validation_root).as_posix(),
                    "sha256": msh_digest,
                },
            }),
            encoding="utf-8",
        )
        native_dir = validation_root / "cgns_to_oneflow/results"
        native_dir.mkdir(parents=True)
        native = native_dir / "lineage.oneflow.grid"
        native.write_bytes(b"native-grid")
        native_digest = hashlib.sha256(native.read_bytes()).hexdigest()
        (native_dir / (native.name + ".geometry-lineage.json")).write_text(
            json.dumps({
                "schema": "oneflow.geometry-artifact-lineage",
                "version": 1,
                "artifact": {
                    "role": "OneFLOW native grid", "file_name": native.name,
                    "format": "GRID", "sha256": native_digest,
                },
                "geometry_lineage": geometry_lineage,
                "parent_artifact": {
                    "role": "CGNS volume mesh", "file_name": cgns.name,
                    "path": cgns.relative_to(validation_root).as_posix(),
                    "sha256": cgns_digest,
                },
            }),
            encoding="utf-8",
        )
        reopened = load_mesh_project(project["directory"])
        self.assertEqual(reopened["geometry_lineage"]["asset_id"], "geometry.sphere.unit")
        self.assertEqual(reopened["geometry_lineage"]["provenance"]["kind"], "generated")
        artifacts = reopened["geometry_lineage"]["downstream_artifacts"]
        by_role = {item["role"]: item for item in artifacts}
        self.assertEqual(
            set(by_role),
            {"CGNS volume mesh", "OneFLOW native grid"},
        )
        self.assertEqual(
            by_role["OneFLOW native grid"]["parent"]["file_name"],
            cgns.name,
        )

        # A replaced CGNS file invalidates both its own lineage and the
        # native-grid record that points to its old checksum.
        cgns.write_bytes(b"replaced-cgns")
        reopened = load_mesh_project(project["directory"])
        self.assertEqual(reopened["geometry_lineage"].get("downstream_artifacts", []), [])

    def test_project_reopen_does_not_trust_mismatched_stl_manifest(self) -> None:
        source = self.home / "untrusted.stl"
        source.write_bytes(b"actual-stl")
        source.with_suffix(".geometry.json").write_text(
            json.dumps({
                "schema": "oneflow.stl-geometry-manifest",
                "asset_id": "geometry.forged",
                "stl_sha256": "0" * 64,
            }),
            encoding="utf-8",
        )
        project = create_mesh_project_from_stl(
            self.repository_root, source, "untrusted_project"
        )
        self.assertNotIn("asset_id", project["geometry_lineage"])
        self.assertEqual(project["geometry_lineage"]["provenance"]["kind"], "unknown")
        self.assertEqual(
            project["geometry_lineage"]["manifest_status"],
            "invalid_or_checksum_mismatch",
        )

    def test_legacy_imported_project_uses_geometry_neutral_mesh_handoff(self) -> None:
        project_root = self.home / "legacy_project"
        input_directory = project_root / "input"
        input_directory.mkdir(parents=True)
        source = self.repository_root / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python/input/model.stl"
        (input_directory / "model.stl").write_bytes(source.read_bytes())
        manifest = {
            "schema": "oneflow.mesh-generation-project",
            "version": 1,
            "name": "legacy_project",
            "template_case_id": None,
            "mesh_case_id": "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            "input_stl": "input/model.stl",
            "output_mesh": "results/legacy_project.msh",
        }
        (project_root / "project.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        loaded = load_mesh_project(project_root)
        expected = "stl_to_gmsh_generic_3d_mesh_handoff"
        self.assertEqual(loaded["mesh_case_id"], expected)
        saved = json.loads((project_root / "project.json").read_text(encoding="utf-8"))
        self.assertEqual(saved["mesh_case_id"], expected)

    def test_import_stl_copies_database_geometry_into_independent_project(self) -> None:
        source = self.repository_root / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python/input/model.stl"
        project = create_mesh_project_from_stl(self.repository_root, source, "imported_spheroid")
        self.assertTrue(project["stl_path"].is_file())
        self.assertNotEqual(project["stl_path"], source)
        self.assertEqual(project["mesh_case_id"], "stl_to_gmsh_generic_3d_mesh_handoff")


if __name__ == "__main__":
    unittest.main()
