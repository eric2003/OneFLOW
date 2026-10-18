import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "tools" / "mesh_layout" / "compile_rectangle.py"
SPEC_PATH = ROOT / "grid" / "minimal_rectangle" / "rectangle.json"

module_spec = importlib.util.spec_from_file_location("compile_rectangle", MODULE_PATH)
compiler = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(compiler)


class CompileRectangleTest(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))

    def test_example_compiles_to_existing_atomic_language(self):
        output = compiler.compile_rectangle(self.spec)
        self.assertIn("Point(1)= {0,0,0}", output)
        self.assertIn("Point(3)= {1,1,0}", output)
        self.assertIn("Line(4)={4,1}", output)
        self.assertIn("Dim(1)={3}", output)
        self.assertIn("Ds(4)={d,-1,-1}", output)
        self.assertIn("Add L2F 1 4 4", output)
        self.assertTrue(output.endswith("Add F2B 1 1 1\n"))

    def test_edge_point_counts_are_independent(self):
        self.spec["edgePointCounts"] = [2, 3, 2, 3]
        output = compiler.compile_rectangle(self.spec)
        for edge, count in enumerate((2, 3, 2, 3), start=1):
            self.assertIn(f"Dim({edge})={{{count}}}", output)

    def test_svg_preview_tracks_point_counts_and_line_ids(self):
        svg = compiler.render_rectangle_svg(self.spec)
        self.assertIn('data-line-id="1"', svg)
        self.assertIn("Line 1: 3 points", svg)
        self.assertIn("Point 4", svg)

    def test_rejects_mismatched_opposite_edge_counts(self):
        self.spec["edgePointCounts"] = [2, 3, 4, 3]
        with self.assertRaisesRegex(ValueError, "opposite edges"):
            compiler.compile_rectangle(self.spec)

    def test_rejects_missing_edge_point_count(self):
        self.spec["edgePointCounts"] = [3, 3, 3]
        with self.assertRaisesRegex(ValueError, "four values"):
            compiler.compile_rectangle(self.spec)

    def test_rejects_invalid_edge_point_count(self):
        self.spec["edgePointCounts"] = [3, 1, 3, 3]
        with self.assertRaisesRegex(ValueError, "integer >= 2"):
            compiler.compile_rectangle(self.spec)

    def test_named_boundary_types_compile_to_legacy_codes(self):
        self.spec["boundaries"] = {
            "line1": "wall",
            "line2": "outflow",
            "line3": "symmetry",
            "line4": "farfield",
        }
        output = compiler.compile_rectangle(self.spec)
        self.assertIn("Boundary(1)={2}", output)
        self.assertIn("Boundary(2)={6}", output)
        self.assertIn("Boundary(3)={3}", output)
        self.assertIn("Boundary(4)={4}", output)

    def test_rejects_unknown_boundary_type(self):
        self.spec["boundaries"] = {"line1": "made_up"}
        with self.assertRaisesRegex(ValueError, "known boundary name"):
            compiler.compile_rectangle(self.spec)


if __name__ == "__main__":
    unittest.main()
