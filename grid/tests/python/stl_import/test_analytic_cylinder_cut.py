"""Tests for the first analytic-cylinder cutting stage."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from grid.stl_import.analytic_cylinder import build_cylinder_lateral_surface


class AnalyticCylinderCutTest(unittest.TestCase):
    def test_two_axial_planes_create_two_open_boundary_loops(self):
        triangles = build_cylinder_lateral_surface(1.0, -2.0, 2.0, -1.5, 1.5, 48, 12)
        self.assertEqual(len(triangles), 48 * 12 * 2)
        for triangle in triangles:
            for point in triangle:
                self.assertAlmostEqual(point[1] ** 2 + point[2] ** 2, 1.0, places=12)
        root = Path(__file__).resolve().parents[4]
        generator = root / "grid" / "stl_import" / "analytic_cylinder.py"
        inspector = root / "grid" / "stl_import" / "inspect.py"
        with tempfile.TemporaryDirectory() as temporary:
            temp = Path(temporary)
            output_stl = temp / "cut" / "cylinder_cut.stl"
            cut_report = temp / "cut" / "cut_report.json"
            inspect_report = temp / "stage_a" / "stage_a_report.json"
            result = subprocess.run([
                sys.executable, str(generator), "--output-stl", str(output_stl),
                "--report", str(cut_report), "--cut-start", "-1.5", "--cut-end", "1.5",
                "--circumferential-divisions", "48", "--axial-divisions", "12"],
                capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
            generated = json.loads(cut_report.read_text(encoding="utf-8"))
            self.assertEqual(generated["source_kind"], "explicit_analytic_definition")
            self.assertFalse(generated["end_caps_generated"])
            self.assertEqual(generated["surface_quality"]["boundary"]["boundary_loop_count"], 2)
            inspected = subprocess.run([
                sys.executable, str(inspector), "--input", str(output_stl),
                "--output", str(inspect_report)],
                capture_output=True, text=True, check=False)
            self.assertEqual(inspected.returncode, 0, msg=inspected.stdout + inspected.stderr)
            quality = json.loads(inspect_report.read_text(encoding="utf-8"))["surface_quality"]
            self.assertEqual(quality["status"], "warning")
            self.assertFalse(quality["closed"])
            self.assertTrue(quality["edge_manifold"])
            self.assertTrue(quality["consistently_oriented"])
            self.assertEqual(quality["boundary_edge_count"], 96)
            self.assertEqual(quality["boundary"]["boundary_component_count"], 2)
            self.assertEqual(quality["boundary"]["boundary_loop_count"], 2)
            self.assertEqual(quality["boundary"]["boundary_component_vertex_counts"], [48, 48])
            self.assertEqual(quality["non_manifold_edge_count"], 0)
            self.assertEqual(quality["orientation_conflict_edge_count"], 0)
            self.assertEqual(quality["connected_component_count"], 1)
            self.assertFalse(quality["volume_mesh_ready"])

    def test_cut_planes_must_be_inside_and_ordered(self):
        for start, end in ((-2.0, 1.0), (-1.0, 2.0), (1.0, -1.0), (0.0, 0.0)):
            with self.subTest(start=start, end=end):
                with self.assertRaises(ValueError):
                    build_cylinder_lateral_surface(1.0, -2.0, 2.0, start, end)


if __name__ == "__main__":
    unittest.main()
