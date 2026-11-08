from __future__ import annotations

import json
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from mesh_cases import load_catalog, validate_case  # noqa: E402
from grid.stl.stl_to_gmsh import convert  # noqa: E402
from tools.oneflow_gui.workflow_plan import build_case_graph_plan  # noqa: E402

EXPECTED_IDS = {
    "stl_to_gmsh_square",
    "stl_to_gmsh_square_2x2",
    "stl_to_gmsh_two_square_cells",
    "stl_to_gmsh_python_square",
    "stl_to_gmsh_python_two_square_cells",
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
        ):
            entry = self.entries[case_id]
            case = ROOT / entry["directory"]
            self.assertNotIn("execution", entry)
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
        ):
            entry = self.entries[case_id]
            self.assertEqual(entry["execution"]["kind"], "python")
            self.assertEqual(validate_case(ROOT, entry), [])

    def test_square_2x2_reuses_square_stl_and_changes_only_mesh_density(self):
        square = ROOT / "grid/stl_to_gmsh/square"
        square_2x2 = ROOT / "grid/stl_to_gmsh/square_2x2"
        self.assertEqual(
            (square / "input/model.stl").read_bytes(),
            (square_2x2 / "input/model.stl").read_bytes(),
        )
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
            ("python/square", 1, 1),
            ("python/two_square_cells", 2, 1),
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
            (ROOT / "grid/stl_to_gmsh/python/square/expected.json").read_text()
        )
        two = json.loads(
            (ROOT / "grid/stl_to_gmsh/python/two_square_cells/expected.json").read_text()
        )
        self.assertEqual(
            square["mesh"],
            {"nodes": 4, "line_elements": 4, "triangle_elements": 0, "quadrilateral_elements": 1},
        )
        self.assertEqual(
            two["mesh"],
            {"nodes": 6, "line_elements": 6, "triangle_elements": 0, "quadrilateral_elements": 2},
        )

    def test_stl_to_cfd_pipeline_declares_ordered_artifact_handoffs(self):
        workflow_path = ROOT / "grid/stl_to_gmsh/develop_test/stl_square_to_cfd.json"
        workflow = json.loads(workflow_path.read_text(encoding="utf-8"))
        plan = build_case_graph_plan(workflow["steps"][0]["cases"], ROOT)
        self.assertEqual(
            [step.case_id for step in plan.steps],
            [
                "stl_to_gmsh_square_2x2",
                "gmsh_to_cgns_square_2x2",
                "cgns_to_oneflow_square_2x2",
                "cfd_solve_square_2x2",
            ],
        )
        gmsh, cgns, native, cfd = plan.steps
        self.assertEqual(gmsh.outputs[0].relative_path, Path("results/square_2x2.msh"))
        self.assertEqual(cgns.inputs[0].producer_id, gmsh.case_id)
        self.assertEqual(cgns.inputs[0].relative_path, Path("results/square_2x2.msh"))
        self.assertEqual(cgns.outputs[0].relative_path, Path("results/square.cgns"))
        self.assertEqual(native.inputs[0].producer_id, cgns.case_id)
        self.assertEqual(native.inputs[0].relative_path, Path("results/square.cgns"))
        self.assertEqual(native.outputs[0].relative_path, Path("results/square.oneflow.grid"))
        self.assertEqual(cfd.inputs[0].producer_id, native.case_id)
        self.assertEqual(cfd.inputs[0].relative_path, Path("results/square.oneflow.grid"))
        self.assertEqual(
            [output.relative_path.as_posix() for output in cfd.outputs],
            ["results/res.dat", "restart/flow.dat", "visual/visual.dat"],
        )
        cgns_grid = (
            ROOT / "grid/stl_to_gmsh/pipeline/cgns_to_oneflow/script/grid.txt"
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
            case_definitions[cfd.case_id]["validation"]["paths"],
            ["results/res.dat", "restart/flow.dat", "visual/visual.dat"],
        )
        cfd_parameters = (
            ROOT / "grid/stl_to_gmsh/pipeline/cfd_uns2dslau2/script/cfd.txt"
        ).read_text(encoding="utf-8")
        self.assertIn("int    maxSteps   =  50;", cfd_parameters)
        self.assertIn("Real global_dt = 0.001;", cfd_parameters)
        self.assertIn("int idump = 0;", cfd_parameters)


if __name__ == "__main__":
    unittest.main()
