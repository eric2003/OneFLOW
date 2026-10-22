"""Tests for the Qt-independent mesh layout editor document service."""

import json
import tempfile
import unittest
from pathlib import Path

from tools.mesh_layout.editor_document import (
    compile_specification,
    load_specification,
    save_specification,
)


ROOT = Path(__file__).resolve().parents[2]


class MeshLayoutEditorDocumentTest(unittest.TestCase):
    def test_rectangle_template_compiles_to_checked_in_layout(self):
        spec_path = ROOT / "grid" / "minimal_rectangle" / "rectangle.json"
        layout_path = ROOT / "grid" / "minimal_rectangle" / "grid" / "minimalRectangleLayout.txt"
        spec = load_specification(spec_path)
        layout, svg = compile_specification(spec)
        self.assertEqual(layout, layout_path.read_text(encoding="utf-8"))
        self.assertIn("Rectangle layout", svg)

    def test_vertical_midpoint_rectangle_template_compiles_to_checked_in_layout(self):
        spec_path = ROOT / "grid" / "minimal_rectangle_vertical_midpoints" / "rectangle_vertical_midpoints.json"
        layout_path = ROOT / "grid" / "minimal_rectangle_vertical_midpoints" / "grid" / "minimalRectangleVerticalMidpointsLayout.txt"
        spec = load_specification(spec_path)
        layout, svg = compile_specification(spec)
        self.assertEqual(layout, layout_path.read_text(encoding="utf-8"))
        self.assertIn("Rectangle with vertical midpoints", svg)
        self.assertEqual(layout.count("Line("), 6)
        self.assertEqual(layout.count("Add L2F"), 6)

    def test_flat_plate_template_compiles_to_checked_in_layout(self):
        spec_path = ROOT / "grid" / "turbulent_flat_plate" / "flat_plate.json"
        layout_path = ROOT / "grid" / "turbulent_flat_plate" / "grid" / "turbulentFlatPlateLayout.txt"
        spec = load_specification(spec_path)
        layout, svg = compile_specification(spec)
        self.assertEqual(layout, layout_path.read_text(encoding="utf-8"))
        self.assertIn("Flat plate boundary-layer layout", svg)

    def test_ten_edge_laminar_flat_plate_compiles_to_checked_in_layout(self):
        spec_path = ROOT / "grid" / "laminar_flat_plate_10_edges" / "flat_plate.json"
        layout_path = ROOT / "grid" / "laminar_flat_plate_10_edges" / "grid" / "laminarFlatPlate10EdgeLayout.txt"
        spec = load_specification(spec_path)
        layout, svg = compile_specification(spec)
        self.assertEqual(layout, layout_path.read_text(encoding="utf-8"))
        self.assertIn("Flat plate boundary-layer layout", svg)
        self.assertIn("41 points", svg)
        self.assertEqual(layout.count("Line("), 10)
        self.assertIn("Dim(3)={41}", layout)
        self.assertIn("Dim(4)={45}", layout)
        self.assertIn("Ds(3)={d,-1,-1}", layout)
        self.assertIn("Ds(4)={tanh,0.0005,-1}", layout)
        self.assertEqual(layout.count("Add L2F"), 12)

    def test_save_specification_writes_valid_ascii_json(self):
        source = ROOT / "grid" / "minimal_rectangle" / "rectangle.json"
        spec = load_specification(source)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "saved" / "rectangle.json"
            save_specification(output, spec)
            raw = output.read_bytes()
            self.assertTrue(raw.endswith(b"\n"))
            self.assertTrue(raw.isascii())
            self.assertEqual(json.loads(raw), spec)

    def test_unsupported_template_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unsupported mesh layout type"):
            compile_specification({"type": "circle"})


if __name__ == "__main__":
    unittest.main()
