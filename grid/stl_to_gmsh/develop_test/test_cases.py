from __future__ import annotations
import json, tempfile, unittest, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
from mesh_cases import load_catalog, validate_case  # noqa: E402
from grid.stl.stl_to_gmsh import convert  # noqa: E402
EXPECTED_IDS = {"stl_to_gmsh_square", "stl_to_gmsh_two_square_cells",
                "stl_to_gmsh_python_square", "stl_to_gmsh_python_two_square_cells"}

class StlToGmshCaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = load_catalog(ROOT / "grid/case_catalog.json")
        cls.entries = {e["id"]: e for e in catalog["cases"]
                       if e.get("group") == "STL to Gmsh"}

    def test_native_and_python_cases_are_registered_separately(self):
        self.assertEqual(set(self.entries), EXPECTED_IDS)
        for case_id in ("stl_to_gmsh_square", "stl_to_gmsh_two_square_cells"):
            entry = self.entries[case_id]
            case = ROOT / entry["directory"]
            self.assertNotIn("execution", entry)
            self.assertEqual(validate_case(ROOT, entry), [])
            self.assertEqual(sorted(p.name for p in (case / "script").glob("*.txt")),
                             ["dimension.txt", "grid.txt", "simu.txt"])
            grid = (case / "script/grid.txt").read_text(encoding="utf-8")
            self.assertIn('string sourceGridType = "stl";', grid)
            self.assertIn('string targetGridType = "gmsh";', grid)
            self.assertIn('string gmshExecutable = "gmsh";', grid)
            self.assertFalse((case / "script/control.txt").exists())
        for case_id in ("stl_to_gmsh_python_square", "stl_to_gmsh_python_two_square_cells"):
            entry = self.entries[case_id]
            self.assertEqual(entry["execution"]["kind"], "python")
            self.assertEqual(validate_case(ROOT, entry), [])

    def test_python_geometry_compiles_square_and_two_cell_layouts(self):
        for relative, nx, ny in (("python/square", 1, 1),
                                 ("python/two_square_cells", 2, 1)):
            case = ROOT / "grid/stl_to_gmsh" / relative
            with tempfile.TemporaryDirectory() as temp:
                geo = Path(temp) / "geometry.geo"
                report = convert(case / "input/model.stl",
                                 case / "script/layout.json", geo, tolerance=1.0e-10)
                self.assertEqual(report["boundary_vertices"], 4)
                layout = json.loads((case / "script/layout.json").read_text(encoding="utf-8"))
                self.assertEqual(layout["transfinite_cells"], {"x": nx, "y": ny})
                text = geo.read_text(encoding="utf-8")
                for name in ("inflow", "outflow", "symmetry", "farfield"):
                    self.assertIn(name, text)
                self.assertIn('Physical Surface("fluid_domain")', text)
            self.assertFalse((case / "results").exists())

    def test_expected_mesh_counts_distinguish_the_two_examples(self):
        square = json.loads((ROOT / "grid/stl_to_gmsh/python/square/expected.json").read_text())
        two = json.loads((ROOT / "grid/stl_to_gmsh/python/two_square_cells/expected.json").read_text())
        self.assertEqual(square["mesh"], {"nodes": 4, "line_elements": 4,
                         "triangle_elements": 0, "quadrilateral_elements": 1})
        self.assertEqual(two["mesh"], {"nodes": 6, "line_elements": 6,
                         "triangle_elements": 0, "quadrilateral_elements": 2})

if __name__ == "__main__":
    unittest.main()
