"""Regression tests for deterministic STL geometry generation."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from grid.stl_to_gmsh.generate_stl_geometry import generate_geometry


class StlGeometryGenerationTest(unittest.TestCase):
    def test_unit_sphere_matches_documented_icosphere_topology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = generate_geometry("sphere", Path(temporary) / "sphere.stl")
            self.assertEqual(result["surface"]["triangle_count"], 1280)
            self.assertEqual(result["surface"]["vertex_count"], 642)
            self.assertTrue(result["surface"]["watertight"])
            self.assertGreater(result["surface"]["signed_volume"], 0.0)
            self.assertEqual(result["patch_triangle_indices"]["wall"], list(range(1280)))
            self.assertEqual(result["units"], "m")
            self.assertEqual(result["stl_sha256"], result["stl_sha256"])
            self.assertTrue(result["manifest"].is_file())

    def test_prolate_spheroid_uses_configured_semi_axes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = generate_geometry(
                "prolate_spheroid",
                Path(temporary) / "prolate.stl",
                semi_axes=(3.0, 0.5, 0.5),
            )
            self.assertTrue(result["surface"]["watertight"])
            self.assertAlmostEqual(result["surface"]["bounds"]["x"][1], 3.0)
            self.assertAlmostEqual(result["surface"]["bounds"]["y"][1], 0.5)
            self.assertAlmostEqual(result["surface"]["bounds"]["z"][1], 0.5)
            self.assertAlmostEqual(result["surface"]["signed_volume"], math.pi, delta=0.15)

    def test_cylinder_is_closed_and_preserves_three_surface_patch_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = generate_geometry(
                "cylinder_internal_flow",
                Path(temporary) / "cylinder.stl",
                radius=0.5,
                length=4.0,
                circumferential_divisions=32,
                axial_divisions=8,
            )
            surface = result["surface"]
            labels = result["patch_triangle_indices"]
            self.assertTrue(surface["watertight"])
            self.assertAlmostEqual(surface["signed_volume"], 32 * 0.5 * 0.5 * math.sin(2 * math.pi / 32) / 2, delta=0.2)
            self.assertEqual(len(labels["inflow"]), 32)
            self.assertEqual(len(labels["outflow"]), 32)
            self.assertEqual(len(labels["wall"]), 2 * 32 * 8)
            self.assertEqual(
                sum(map(len, labels.values())),
                surface["triangle_count"],
            )
            manifest = json.loads(result["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["physical_groups"]["surface"], ["wall", "inflow", "outflow"])
            self.assertIn("not solid walls", manifest["notes"])

    def test_invalid_geometry_parameters_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaises(ValueError):
                generate_geometry("sphere", Path(temporary) / "bad.stl", radius=0.0)
            with self.assertRaises(ValueError):
                generate_geometry(
                    "cylinder_internal_flow",
                    Path(temporary) / "bad.stl",
                    circumferential_divisions=6,
                )


if __name__ == "__main__":
    unittest.main()
