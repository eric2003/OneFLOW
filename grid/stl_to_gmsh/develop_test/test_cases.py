from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from mesh_cases import load_catalog, validate_case  # noqa: E402
from grid.stl.stl_to_gmsh import convert, inspect_msh4  # noqa: E402
from tools.oneflow_gui.workflow_plan import build_case_graph_plan  # noqa: E402

EXPECTED_IDS = {
    "stl_to_gmsh_square",
    "stl_to_gmsh_square_2x2",
    "stl_to_gmsh_two_square_cells",
    "stl_to_gmsh_circle_external_flow",
    "stl_to_gmsh_python_square",
    "stl_to_gmsh_python_two_square_cells",
    "stl_to_gmsh_python_circle_external_flow",
    "stl_to_gmsh_python_channel_internal_flow",
}


class StlToGmshCaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = load_catalog(ROOT / "grid/case_catalog.json")
        cls.entries = {
            entry["id"]: entry
            for entry in catalog["cases"]
            if entry.get("group") == "STL to Gmsh"
        }

    def test_native_and_python_cases_are_registered_separately(self):
        self.assertEqual(set(self.entries), EXPECTED_IDS)
        for case_id in (
            "stl_to_gmsh_square",
            "stl_to_gmsh_square_2x2",
            "stl_to_gmsh_two_square_cells",
            "stl_to_gmsh_circle_external_flow",
        ):
            entry = self.entries[case_id]
            case = ROOT / entry["directory"]
            self.assertNotIn("execution", entry)
            self.assertTrue(entry["directory"].endswith("/cpp"))
            self.assertEqual(validate_case(ROOT, entry), [])
            self.assertEqual(
                sorted(path.name for path in (case / "script").glob("*.txt")),
                ["dimension.txt", "grid.txt", "simu.txt"],
            )
            grid = (case / "script/grid.txt").read_text(encoding="utf-8")
            self.assertIn('string sourceGridType = "stl";', grid)
            self.assertIn('string targetGridType = "gmsh";', grid)
            self.assertIn('string gmshExecutable = "gmsh";', grid)
            self.assertFalse((case / "script/control.txt").exists())
        for case_id in (
            "stl_to_gmsh_python_square",
            "stl_to_gmsh_python_two_square_cells",
            "stl_to_gmsh_python_circle_external_flow",
            "stl_to_gmsh_python_channel_internal_flow",
        ):
            entry = self.entries[case_id]
            self.assertEqual(entry["execution"]["kind"], "python")
            self.assertTrue(entry["directory"].endswith("/python"))
            self.assertEqual(entry["execution"]["script"], "grid/stl_to_gmsh/run_python_case.py")
            self.assertEqual(validate_case(ROOT, entry), [])

    def test_native_circle_external_flow_uses_same_stl_and_domain_as_python_reference(self):
        native = ROOT / "grid/stl_to_gmsh/circle_external_flow/cpp"
        reference = ROOT / "grid/stl_to_gmsh/circle_external_flow/python"
        self.assertEqual(
            (native / "input/model.stl").read_bytes(),
            (reference / "input/model.stl").read_bytes(),
        )
        native_layout = json.loads(
            (native / "script/layout.json").read_text(encoding="utf-8")
        )
        python_layout = json.loads(
            (reference / "script/layout.json").read_text(encoding="utf-8")
        )
        self.assertEqual(native_layout, python_layout)
        self.assertEqual(native_layout["surface_name"], "fluid_domain")
        self.assertEqual(native_layout["body_boundary_name"], "wall")
        self.assertEqual(
            native_layout["domain"],
            {
                "xmin": -3, "xmax": 6, "ymin": -3, "ymax": 3,
                "boundaries": {
                    "left": "inflow", "right": "outflow",
                    "bottom": "farfield", "top": "farfield",
                },
            },
        )
        self.assertEqual(
            native_layout["quality"],
            {"minimum_triangle_quality": 0.01, "minimum_angle_degrees": 1.0},
        )
        self.assertEqual(native_layout["domain"]["boundaries"], {
            "left": "inflow", "right": "outflow",
            "bottom": "farfield", "top": "farfield",
        })
        grid = (native / "script/grid.txt").read_text(encoding="utf-8")
        self.assertIn('string targetGridFileName = "results/circle_external_flow.msh";', grid)
        self.assertFalse((native / "results").exists())

    def test_square_2x2_reuses_square_stl_and_changes_only_mesh_density(self):
        square = ROOT / "grid/stl_to_gmsh/square/cpp"
        square_2x2 = ROOT / "grid/stl_to_gmsh/square_2x2/cpp"
        self.assertEqual(
            (square / "input/model.stl").read_bytes(),
            (square_2x2 / "input/model.stl").read_bytes(),
        )
        self.assertTrue((square_2x2.parent / "README.md").is_file())
        layout = json.loads(
            (square_2x2 / "script/layout.json").read_text(encoding="utf-8")
        )
        self.assertEqual(layout["transfinite_cells"], {"x": 2, "y": 2})
        grid = (square_2x2 / "script/grid.txt").read_text(encoding="utf-8")
        self.assertIn('string sourceGridFileName = "input/model.stl";', grid)
        self.assertIn(
            'string targetGridFileName = "results/square_2x2.msh";',
            grid,
        )
        self.assertFalse((square_2x2 / "results").exists())

    def test_python_geometry_compiles_square_and_two_cell_layouts(self):
        for relative, nx, ny in (
            ("square/python", 1, 1),
            ("two_square_cells/python", 2, 1),
        ):
            case = ROOT / "grid/stl_to_gmsh" / relative
            with tempfile.TemporaryDirectory() as temp:
                geo = Path(temp) / "geometry.geo"
                report = convert(
                    case / "input/model.stl",
                    case / "script/layout.json",
                    geo,
                    tolerance=1.0e-10,
                )
                self.assertEqual(report["boundary_vertices"], 4)
                layout = json.loads(
                    (case / "script/layout.json").read_text(encoding="utf-8")
                )
                self.assertEqual(layout["transfinite_cells"], {"x": nx, "y": ny})
                text = geo.read_text(encoding="utf-8")
                for name in ("inflow", "outflow", "symmetry", "farfield"):
                    self.assertIn(name, text)
                self.assertIn('Physical Surface("fluid_domain")', text)
            self.assertFalse((case / "results").exists())

    def test_expected_mesh_counts_distinguish_the_two_examples(self):
        square = json.loads(
            (ROOT / "grid/stl_to_gmsh/square/python/expected.json").read_text()
        )
        two = json.loads(
            (ROOT / "grid/stl_to_gmsh/two_square_cells/python/expected.json").read_text()
        )
        self.assertEqual(
            square["mesh"],
            {"nodes": 4, "line_elements": 4, "triangle_elements": 0, "quadrilateral_elements": 1},
        )
        self.assertEqual(
            two["mesh"],
            {"nodes": 6, "line_elements": 6, "triangle_elements": 0, "quadrilateral_elements": 2},
        )

    def test_msh_inspector_measures_triangle_quality_and_domain_bounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            msh = Path(temporary) / "one_triangle.msh"
            msh.write_text(
                "$MeshFormat\\n4.1 0 8\\n$EndMeshFormat\\n"
                '$PhysicalNames\n1\n2 1 "fluid_domain"\n$EndPhysicalNames\n'
                "$Nodes\\n1 3 1 3\\n2 1 0 3\\n1 2 3\\n"
                "0 0 0\\n1 0 0\\n0 1 0\\n$EndNodes\\n"
                "$Elements\\n1 1 1 1\\n2 1 2 1\\n1 1 2 3\\n$EndElements\\n",
                encoding="utf-8",
            )
            summary = inspect_msh4(msh)
        quality = summary["triangle_quality"]
        self.assertEqual(summary["triangle_elements"], 1)
        self.assertEqual(quality["degenerate_count"], 0)
        self.assertAlmostEqual(quality["minimum_area"], 0.5)
        self.assertAlmostEqual(quality["minimum_quality"], 3.0 ** 0.5 / 2.0)
        self.assertAlmostEqual(quality["minimum_angle_degrees"], 45.0)
        self.assertEqual(summary["coordinate_bounds"]["minimum"], [0.0, 0.0, 0.0])
        self.assertEqual(summary["coordinate_bounds"]["maximum"], [1.0, 1.0, 0.0])

    def test_python_circle_external_flow_generates_valid_gmsh_mesh(self):
        gmsh = shutil.which("gmsh")
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/circle_external_flow/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "circle_external_flow"
            shutil.copytree(source, staged)
            command = [
                sys.executable,
                str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory",
                str(staged),
                "--gmsh",
                gmsh,
            ]
            result = subprocess.run(
                command, cwd=ROOT, text=True, capture_output=True, check=False
            )
            self.assertEqual(
                result.returncode, 0,
                msg=result.stdout + "\\n" + result.stderr,
            )
            mesh = inspect_msh4(staged / "results/model.msh")
            self.assertGreater(mesh["nodes"], 20)
            self.assertGreater(mesh["line_elements"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            self.assertEqual(mesh["quadrilateral_elements"], 0)
            quality = mesh["triangle_quality"]
            self.assertEqual(quality["degenerate_count"], 0)
            self.assertGreaterEqual(quality["minimum_quality"], 0.01)
            self.assertGreaterEqual(quality["minimum_angle_degrees"], 1.0)
            self.assertEqual(mesh["coordinate_bounds"]["minimum"][:2], [-3.0, -3.0])
            self.assertEqual(mesh["coordinate_bounds"]["maximum"][:2], [6.0, 3.0])
            groups = {
                (group["dimension"], group["name"])
                for group in mesh["physical_groups"]
            }
            self.assertEqual(groups, {
                (1, "wall"), (1, "inflow"), (1, "outflow"),
                (1, "farfield"), (2, "fluid_domain"),
            })
            report = json.loads(
                (staged / "results/report.json").read_text(encoding="utf-8")
            )
            self.assertTrue(report["validated"])
        self.assertFalse((source / "results").exists())

    def test_external_and_internal_flow_cases_build_distinct_fluid_domains(self):
        cases_root = ROOT / "grid/stl_to_gmsh"
        external = cases_root / "circle_external_flow/python"
        internal = cases_root / "channel_internal_flow/python"
        with tempfile.TemporaryDirectory() as temp:
            external_geo = Path(temp) / "external.geo"
            external_report = convert(
                external / "input/model.stl",
                external / "script/layout.json",
                external_geo,
                tolerance=1.0e-10,
            )
            external_text = external_geo.read_text(encoding="utf-8")
            self.assertEqual(external_report["body_boundary_edges"], 32)
            self.assertIn("Plane Surface(1) = {1, 2};", external_text)
            self.assertIn("lc = 0.22;", external_text)
            self.assertIn("lc_outer =", external_text)
            self.assertIn("Point(33) = {-3, -3, 0, lc_outer};", external_text)
            self.assertIn('Physical Curve("wall")', external_text)
            self.assertEqual(external_text.count('Physical Curve("farfield")'), 1)
            for name in ("inflow", "outflow", "farfield", "wall"):
                self.assertIn(f'Physical Curve("{name}")', external_text)
            self.assertIn('Physical Surface("fluid_domain")', external_text)

            internal_geo = Path(temp) / "internal.geo"
            internal_report = convert(
                internal / "input/model.stl",
                internal / "script/layout.json",
                internal_geo,
                tolerance=1.0e-10,
            )
            internal_text = internal_geo.read_text(encoding="utf-8")
            self.assertEqual(internal_report["boundary_edges"], 4)
            self.assertIn("Plane Surface(1) = {1};", internal_text)
            self.assertNotIn("Plane Surface(1) = {1, 2};", internal_text)
            for name in ("Inlet", "Outlet", "UpperWall", "LowerWall"):
                self.assertIn(f'Physical Curve("{name}")', internal_text)
            self.assertIn("Transfinite Curve", internal_text)
            self.assertIn("Transfinite Surface {1};", internal_text)
            self.assertIn("Recombine Surface {1};", internal_text)

        for case in (external, internal):
            self.assertFalse((case / "results").exists())

    def test_stl_to_cfd_pipeline_declares_ordered_artifact_handoffs(self):
        workflow_path = ROOT / "grid/stl_to_gmsh/develop_test/stl_square_to_cfd.json"
        workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
        plan = build_case_graph_plan(workflow["steps"][0]["cases"], ROOT)
        self.assertEqual(
            [step.case_id for step in plan.steps],
            [
                "stl_to_gmsh_square_2x2",
                "gmsh_to_cgns_square_2x2",
                "gmsh_to_oneflow_square_2x2",
                "cfd_solve_square_2x2",
            ],
        )
        gmsh, cgns, native, cfd = plan.steps
        self.assertEqual(gmsh.outputs[0].relative_path, Path("results/square_2x2.msh"))
        self.assertEqual(
            gmsh.outputs[1].relative_path,
            Path("results/square_2x2.stl-to-gmsh-report.json"),
        )
        self.assertEqual(cgns.inputs[0].producer_id, gmsh.case_id)
        self.assertEqual(cgns.inputs[0].relative_path, Path("results/square_2x2.msh"))
        self.assertEqual(cgns.outputs[0].relative_path, Path("results/square.cgns"))
        self.assertEqual(
            [output.relative_path.as_posix() for output in cgns.outputs],
            [
                "results/square.cgns",
                "results/square.gmsh-read-report.json",
                "results/square.gmsh-report.md",
            ],
        )
        self.assertEqual(native.inputs[0].producer_id, cgns.case_id)
        self.assertEqual(native.inputs[0].relative_path, Path("results/square.cgns"))
        self.assertEqual(native.outputs[0].relative_path, Path("results/square.oneflow.grid"))
        self.assertEqual(native.outputs[1].relative_path, Path("results/bcRegionMap.txt"))
        self.assertEqual(native.outputs[2].relative_path, Path("results/square.info.json"))
        self.assertEqual(cfd.inputs[0].producer_id, native.case_id)
        self.assertEqual(cfd.inputs[0].relative_path, Path("results/square.oneflow.grid"))
        self.assertEqual(
            [output.relative_path.as_posix() for output in cfd.outputs],
            ["results/res.dat", "restart/flow.dat", "visual/visual.dat"],
        )
        cgns_grid = (
            ROOT / "grid/stl_to_gmsh/pipeline/gmsh_to_oneflow/script/grid.txt"
        ).read_text(encoding="utf-8")
        self.assertIn('string sourceGridType = "cgns";', cgns_grid)
        self.assertIn('string targetGridType = "oneflow";', cgns_grid)
        case_definitions = {
            case["id"]: case for case in workflow["steps"][0]["cases"]
        }
        self.assertEqual(case_definitions[gmsh.case_id]["validation"]["kind"], "gmsh_msh4")
        self.assertEqual(
            case_definitions[gmsh.case_id]["validation"]["expected_counts"],
            {
                "elements": 12,
                "nodes": 9,
                "line_elements": 8,
                "triangle_elements": 0,
                "quadrilateral_elements": 4,
            },
        )
        self.assertEqual(
            case_definitions[gmsh.case_id]["validation"]["expected_element_types"],
            {"1": 8, "3": 4},
        )
        for case_id in (cgns.case_id, native.case_id, cfd.case_id):
            self.assertEqual(
                case_definitions[case_id]["validation"]["kind"],
                "non_empty_files",
            )
        self.assertEqual(
            case_definitions[cgns.case_id]["validation"]["paths"],
            [
                "results/square.cgns",
                "results/square.gmsh-read-report.json",
                "results/square.gmsh-report.md",
            ],
        )
        self.assertEqual(
            case_definitions[native.case_id]["validation"]["paths"],
            [
                "results/square.oneflow.grid",
                "results/bcRegionMap.txt",
                "results/square.info.json",
            ],
        )
        self.assertEqual(
            case_definitions[cfd.case_id]["validation"]["paths"],
            ["results/res.dat", "restart/flow.dat", "visual/visual.dat"],
        )
        cfd_parameters = (
            ROOT / "grid/stl_to_gmsh/pipeline/cfd_uns2dslau2/script/cfd.txt"
        ).read_text(encoding="utf-8")
        self.assertIn("int    maxSteps   =  50;", cfd_parameters)
        self.assertIn("Real global_dt = 0.001;", cfd_parameters)
        self.assertIn("int idump = 0;", cfd_parameters)
        self.assertIn('string gridFileName = "grid/square.oneflow.grid";', cfd_parameters)


if __name__ == "__main__":
    unittest.main()
