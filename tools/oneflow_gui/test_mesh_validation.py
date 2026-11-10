"""Mesh-to-CFD validation Case catalog and graph tests."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from tools.oneflow_gui.mesh_validation import (
    build_mesh_validation_cases, cfd_template_choices, mesh_case_choices,
    mesh_validation_workspace, reset_mesh_validation_workspace,
    _graph_step_records, _mesh_case_command,
)
from tools.oneflow_gui.workflow_graph_runner import (
    CaseExecutionResult, WorkflowGraphRunResult,
)


class MeshValidationCatalogTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "grid/stl_to_gmsh/square/single_quad/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/square/two_triangles/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/square/two_triangles/python").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/circle_external_flow/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/cube_single_hex/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/square/single_quad/python").mkdir(parents=True)
        (self.root / "test/cfd/plateuns2dslau2/script").mkdir(parents=True)
        (self.root / "grid/case_catalog.json").write_text(json.dumps({
            "schema": "oneflow.case-catalog",
            "cases": [
                {"id": "stl_to_gmsh_square", "title": "One Square QUAD_4 Cell",
                 "directory": "grid/stl_to_gmsh/square/single_quad/cpp", "group": "STL to Gmsh"},
                {"id": "stl_to_gmsh_square_two_triangles",
                 "title": "Square with Two TRI_3 Cells (Native C++)",
                 "directory": "grid/stl_to_gmsh/square/two_triangles/cpp",
                 "group": "STL to Gmsh"},
                {"id": "stl_to_gmsh_circle_external_flow",
                 "title": "Circle External Flow (Native C++)",
                 "directory": "grid/stl_to_gmsh/circle_external_flow/cpp",
                 "group": "STL to Gmsh"},
                {"id": "native_cube", "title": "Cube Single Hex (Native C++)",
                 "directory": "grid/stl_to_gmsh/cube_single_hex/cpp",
                 "group": "STL to Gmsh"},
                {"id": "python_cube", "title": "Cube Single Hex (Python Reference)",
                 "directory": "grid/stl_to_gmsh/cube_single_hex/python",
                 "group": "STL to Gmsh",
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
                {"id": "python_square", "title": "One Square QUAD_4 Cell (Python Reference)",
                 "directory": "grid/stl_to_gmsh/square/single_quad/python",
                 "group": "STL to Gmsh",
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
                {"id": "stl_to_gmsh_python_square_two_triangles",
                 "title": "Square with Two TRI_3 Cells (Python Reference)",
                 "directory": "grid/stl_to_gmsh/square/two_triangles/python",
                 "group": "STL to Gmsh",
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
            ],
        }), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/square/single_quad/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "square.msh"}}), encoding="utf-8"
        )
        (self.root / "grid/stl_to_gmsh/circle_external_flow/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/circle_external_flow.msh";\n',
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/cube_single_hex/cpp/script/grid.txt").write_text(
            'string sourceGridFileName = "input/model.stl";\n'
            'string targetGridFileName = "results/cube_single_hex.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/cube_single_hex/cpp/script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/square/single_quad/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/square.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/square/two_triangles/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/square_two_triangles.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/square/two_triangles/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "square_two_triangles.msh"}}),
            encoding="utf-8")
        (self.root / "test/cfd/plateuns2dslau2/script/cfd.txt").write_text(
            'string resFile = "results/res.dat";\n'
            'string flowFile = "restart/flow.dat";\n'
            'string visualFile = "visual/visual.dat";\n', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_mesh_selector_contains_native_and_python_implementations(self):
        self.assertEqual(mesh_case_choices(self.root), (
            ("stl_to_gmsh_square", "One Square QUAD_4 Cell"),
            ("stl_to_gmsh_square_two_triangles", "Square with Two TRI_3 Cells (Native C++)"),
            ("stl_to_gmsh_circle_external_flow", "Circle External Flow (Native C++)"),
            ("native_cube", "Cube Single Hex (Native C++)"),
            ("python_cube", "Cube Single Hex (Python Reference)"),
            ("python_square", "One Square QUAD_4 Cell (Python Reference)"),
            ("stl_to_gmsh_python_square_two_triangles", "Square with Two TRI_3 Cells (Python Reference)"),
        ))
    def test_mesh_producer_command_dispatches_by_implementation(self):
        python_entry = {
            "execution": {
                "kind": "python",
                "script": "grid/stl_to_gmsh/run_python_case.py",
            }
        }
        command = _mesh_case_command(
            self.root, python_entry, Path("OneFLOW"), Path("gmsh"),
            self.root / ".oneflow/staged/python_case",
        )
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(command[1], str((self.root / "grid/stl_to_gmsh/run_python_case.py").resolve()))
        self.assertEqual(command[2:], [
            "--case-directory", str(self.root / ".oneflow/staged/python_case"),
            "--gmsh", "gmsh",
        ])
        native = _mesh_case_command(
            self.root, {}, Path("OneFLOW"), Path("gmsh"),
            self.root / ".oneflow/staged/native_case",
        )
        self.assertEqual(native, [
            "OneFLOW", "r", str(self.root / ".oneflow/staged/native_case")
        ])

    def test_python_mesh_producer_script_cannot_escape_repository(self):
        with self.assertRaisesRegex(ValueError, "inside the repository"):
            _mesh_case_command(
                self.root,
                {"execution": {"kind": "python", "script": "../outside.py"}},
                Path("OneFLOW"), Path("gmsh"), self.root / "staged",
            )

    def test_circle_external_flow_is_available_to_mesh_to_cfd_validation(self):
        cases = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_circle_external_flow", "laminar_flat_plate"
        )
        self.assertEqual(len(cases), 4)
        self.assertEqual(
            cases[0]["directory"], "grid/stl_to_gmsh/circle_external_flow/cpp"
        )
        self.assertEqual(
            cases[0]["outputs"][0]["path"],
            "results/circle_external_flow.msh",
        )
        self.assertEqual(cases[1]["id"], "gmsh_to_cgns")
        self.assertEqual(
            cases[1]["outputs"][0]["path"],
            "results/circle_external_flow.cgns",
        )
        self.assertEqual(cases[2]["id"], "cgns_to_oneflow")
        self.assertEqual(
            cases[2]["outputs"][0]["path"],
            "results/circle_external_flow.oneflow.grid",
        )
        self.assertEqual(cases[2]["outputs"][1], {
            "id": "bc_region_map", "path": "results/circle_external_flow_bcRegionMap.txt",
            "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile",
        })
        self.assertEqual(cases[2]["outputs"][2], {
            "id": "grid_info", "path": "results/circle_external_flow.info.json",
            "parameter_file": "script/grid.txt", "setting": "gridInfoFileName",
        })
        self.assertEqual(cases[3]["id"], "oneflow_cfd_solve")
        self.assertEqual(
            cases[3]["inputs"],
            [
                {"artifact": "oneflow_grid", "setting": "gridFileName",
                 "parameter_file": "script/cfd.txt"},
                {"artifact": "bc_region_map", "setting": "bcRegionMapFile",
                 "parameter_file": "script/cfd.txt"},
            ],
        )

    def test_square_two_triangles_cases_are_available_to_mesh_to_cfd_validation(self):
        native = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_square_two_triangles", "laminar_flat_plate"
        )
        python = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_python_square_two_triangles", "laminar_flat_plate"
        )
        self.assertEqual(len(native), 4)
        self.assertEqual(len(python), 4)
        self.assertEqual(
            native[0]["directory"], "grid/stl_to_gmsh/square/two_triangles/cpp"
        )
        self.assertEqual(
            native[0]["outputs"][0]["path"], "results/square_two_triangles.msh"
        )
        self.assertEqual(
            python[0]["directory"], "grid/stl_to_gmsh/square/two_triangles/python"
        )
        self.assertEqual(
            python[0]["outputs"][0]["path"], "results/square_two_triangles.msh"
        )
        self.assertEqual(native[1]["id"], "gmsh_to_cgns")
        self.assertEqual(python[1]["id"], "gmsh_to_cgns")

    def test_python_mesh_case_uses_the_same_artifact_handoff_contract(self):
        cases = build_mesh_validation_cases(
            self.root, "python_square", "laminar_flat_plate"
        )
        self.assertEqual(len(cases), 4)
        self.assertEqual(cases[0]["directory"], "grid/stl_to_gmsh/square/single_quad/python")
        self.assertEqual(cases[0]["outputs"], [
            {"id": "gmsh_mesh", "path": "results/square.msh"}
        ])
        self.assertEqual(cases[0]["validation"], {
            "kind": "gmsh_msh4", "path": "results/square.msh"
        })
        self.assertEqual(cases[1]["inputs"][0]["artifact"], "gmsh_mesh")
        self.assertEqual(cases[1]["inputs"][0]["setting"], "sourceGridFileName")
        self.assertEqual(cases[1]["directory"], "grid/conversion/gmsh_to_cgns")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")
        self.assertEqual(cases[3]["inputs"][1]["setting"], "bcRegionMapFile")


    def test_native_cube_single_hex_runs_the_complete_cfd_smoke_chain(self):
        cases = build_mesh_validation_cases(self.root, "native_cube", "laminar_flat_plate")
        self.assertEqual([item["id"] for item in cases],
                         ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
                          "oneflow_cfd_solve"])
        self.assertEqual(cases[0]["directory"], "grid/stl_to_gmsh/cube_single_hex/cpp")
        self.assertEqual(cases[0]["runtime_directory"], "stl_to_gmsh")
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["inputs"][0]["artifact"], "oneflow_grid")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")
        self.assertEqual(cases[0]["outputs"][0]["path"], "results/cube_single_hex.msh")
        self.assertEqual(cases[2]["outputs"], [
            {"id": "oneflow_grid", "path": "results/cube_single_hex.oneflow.grid",
             "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
            {"id": "bc_region_map", "path": "results/cube_single_hex_bcRegionMap.txt",
             "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile"},
            {"id": "grid_info", "path": "results/cube_single_hex.info.json",
             "parameter_file": "script/grid.txt", "setting": "gridInfoFileName"},
        ])

    def test_python_cube_single_hex_runs_the_complete_cfd_smoke_chain(self):
        cube_dir = self.root / "grid/stl_to_gmsh/cube_single_hex/python"
        (cube_dir / "script").mkdir(parents=True, exist_ok=True)
        (cube_dir / "script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8"
        )
        (cube_dir / "expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "cube_single_hex.msh"}}),
            encoding="utf-8",
        )
        cases = build_mesh_validation_cases(
            self.root, "python_cube", "laminar_flat_plate"
        )
        self.assertEqual([item["id"] for item in cases],
                         ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
                          "oneflow_cfd_solve"])
        self.assertEqual(cases[0]["validation"], {
            "kind": "gmsh_msh4", "path": "results/cube_single_hex.msh"
        })
        self.assertEqual(cases[2]["outputs"][1]["path"],
                         "results/cube_single_hex_bcRegionMap.txt")
        self.assertEqual(cases[2]["outputs"][2]["path"],
                         "results/cube_single_hex.info.json")
        self.assertIn("oneflow_cfd_solve", [item["id"] for item in cases])

    def test_first_cfd_template_is_the_existing_laminar_flat_plate_case(self):
        self.assertEqual(cfd_template_choices(), (("laminar_flat_plate", "Laminar Flat Plate"),))

    def test_selected_case_builds_conversion_chain_and_real_cfd_consumer(self):
        cases = build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "laminar_flat_plate")
        self.assertEqual(len(cases), 4)
        self.assertEqual(cases[0]["id"], "stl_to_gmsh")
        self.assertEqual(cases[0]["outputs"][0]["id"], "gmsh_mesh")
        self.assertEqual(cases[0]["outputs"][0]["path"], "results/square.msh")
        self.assertEqual(
            [item["runtime_directory"] for item in cases],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(cases[1]["id"], "gmsh_to_cgns")
        self.assertEqual(cases[1]["depends_on"], ["stl_to_gmsh"])
        self.assertEqual(cases[1]["directory"], "grid/conversion/gmsh_to_cgns")
        self.assertEqual(cases[1]["outputs"][0], {
            "id": "cgns_mesh", "path": "results/square.cgns",
            "parameter_file": "script/grid.txt", "setting": "targetGridFileName",
        })
        self.assertEqual(cases[2]["id"], "cgns_to_oneflow")
        self.assertEqual(cases[2]["depends_on"], ["gmsh_to_cgns"])
        self.assertEqual(cases[2]["directory"], "grid/conversion/cgns_to_oneflow")
        self.assertEqual(cases[2]["outputs"], [
            {"id": "oneflow_grid", "path": "results/square.oneflow.grid",
             "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
            {"id": "bc_region_map", "path": "results/square_bcRegionMap.txt",
             "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile"},
            {"id": "grid_info", "path": "results/square.info.json",
             "parameter_file": "script/grid.txt", "setting": "gridInfoFileName"},
        ])
        self.assertEqual(cases[3]["id"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["depends_on"], ["cgns_to_oneflow"])
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["inputs"], [
            {"artifact": "oneflow_grid", "setting": "gridFileName",
             "parameter_file": "script/cfd.txt"},
            {"artifact": "bc_region_map", "setting": "bcRegionMapFile",
             "parameter_file": "script/cfd.txt"},
        ])

    def test_report_records_final_graph_status_and_blocked_steps(self):
        result = WorkflowGraphRunResult("failed", (
            CaseExecutionResult("stl_to_gmsh", "succeeded", "OneFLOW exited successfully."),
            CaseExecutionResult(
                "gmsh_to_cgns", "failed",
                "Declared output artifacts were not produced: results/mesh.cgns",
            ),
            CaseExecutionResult(
                "cgns_to_oneflow", "skipped", "Blocked by unsuccessful dependencies: gmsh_to_cgns",
            ),
            CaseExecutionResult(
                "oneflow_cfd_solve", "skipped", "Blocked by unsuccessful dependencies: cgns_to_oneflow",
            ),
        ))
        self.assertEqual(_graph_step_records(result), [
            {"case_id": "stl_to_gmsh", "status": "succeeded",
             "message": "OneFLOW exited successfully."},
            {"case_id": "gmsh_to_cgns", "status": "failed",
             "message": "Declared output artifacts were not produced: results/mesh.cgns"},
            {"case_id": "cgns_to_oneflow", "status": "skipped",
             "message": "Blocked by unsuccessful dependencies: gmsh_to_cgns"},
            {"case_id": "oneflow_cfd_solve", "status": "skipped",
             "message": "Blocked by unsuccessful dependencies: cgns_to_oneflow"},
        ])

    def test_unknown_mesh_case_and_template_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown mesh-generation Case"):
            build_mesh_validation_cases(self.root, "missing", "laminar_flat_plate")
        with self.assertRaisesRegex(ValueError, "Unknown CFD template"):
            build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "missing")



    def test_workspace_uses_mesh_case_relative_directory(self):
        paths = type(
            "Paths", (),
            {"mesh_validation": self.root / ".oneflow" / "mesh-validation"},
        )()
        with mock.patch(
            "tools.oneflow_gui.mesh_validation.UserDataPaths.discover",
            return_value=paths,
        ):
            workspace = mesh_validation_workspace(self.root, "stl_to_gmsh_square")
        self.assertEqual(
            workspace,
            self.root / ".oneflow" / "mesh-validation" / "stl_to_gmsh" / "square" / "single_quad" / "cpp",
        )

    def test_workspace_rejects_cases_outside_grid(self):
        catalog = json.loads(
            (self.root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        catalog["cases"][0]["directory"] = "../outside"
        (self.root / "grid/case_catalog.json").write_text(
            json.dumps(catalog), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "below grid"):
            mesh_validation_workspace(self.root, "stl_to_gmsh_square")

    def test_reset_clears_only_the_selected_case_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "mesh-validation/stl_to_gmsh/square"
            sibling = root / "mesh-validation/stl_to_gmsh/square_2x2"
            selected.mkdir(parents=True)
            sibling.mkdir(parents=True)
            (selected / "stale.txt").write_text("old", encoding="utf-8")
            (sibling / "keep.txt").write_text("keep", encoding="utf-8")
            reset_mesh_validation_workspace(selected)
            self.assertFalse((selected / "stale.txt").exists())
            self.assertTrue((sibling / "keep.txt").is_file())


if __name__ == "__main__":
    unittest.main()
