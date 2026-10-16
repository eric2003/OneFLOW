import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "tools" / "mesh_layout" / "compile_flat_plate.py"
SPEC_PATH = ROOT / "grid" / "laminar_flat_plate" / "flat_plate.json"

module_spec = importlib.util.spec_from_file_location("compile_flat_plate", MODULE_PATH)
compiler = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(compiler)


class CompileFlatPlateTest(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads(SPEC_PATH.read_text(encoding="utf-8"))

    def test_validation_case_matches_reference_dimensions(self):
        compiler.validate_flat_plate_spec(self.spec)
        self.assertEqual(self.spec["upstreamPointCount"], 33)
        self.assertEqual(self.spec["platePointCount"], 41)
        self.assertEqual(self.spec["normalPointCount"], 85)

    def test_compiles_two_block_layout(self):
        output = compiler.compile_flat_plate(self.spec)
        self.assertIn("Point(1)= {-2,0,0}", output)
        self.assertIn("Point(3)= {2,0,0}", output)
        self.assertIn("Dim(1)={33}", output)
        self.assertIn("Dim(2)={41}", output)
        self.assertIn("Dim(3)={85}", output)
        self.assertIn("Boundary(1)={3}", output)
        self.assertIn("Boundary(2)={2}", output)
        self.assertIn("Boundary(4)={-1}", output)
        self.assertIn("Boundary(5)={6}", output)
        self.assertIn("Add F2B 2 1 2", output)

    def test_compiled_layout_has_expected_mesh_size(self):
        self.assertEqual(
            (self.spec["upstreamPointCount"] + self.spec["platePointCount"] - 1)
            * self.spec["normalPointCount"],
            73 * 85,
        )
        self.assertEqual(
            (self.spec["upstreamPointCount"] - 1)
            + (self.spec["platePointCount"] - 1),
            72,
        )

    def test_boundary_names_map_to_oneflow_codes(self):
        boundaries = self.spec["boundaries"]
        self.assertEqual(compiler._boundary_type(boundaries, "plate"), 2)
        self.assertEqual(compiler._boundary_type(boundaries, "left"), 4)
        self.assertEqual(compiler._boundary_type(boundaries, "right"), 6)
        self.assertEqual(compiler._boundary_type(boundaries, "interface"), -1)

    def test_svg_contains_two_blocks_and_boundary_roles(self):
        svg = compiler.render_flat_plate_svg(self.spec)
        self.assertIn("block interface", svg)
        self.assertIn("farfield", svg)
        self.assertIn("outflow", svg)
        self.assertIn("wall", svg)

    def test_rejects_invalid_normal_point_count(self):
        self.spec["normalPointCount"] = 1
        with self.assertRaisesRegex(ValueError, "normalPointCount"):
            compiler.compile_flat_plate(self.spec)


if __name__ == "__main__":
    unittest.main()
