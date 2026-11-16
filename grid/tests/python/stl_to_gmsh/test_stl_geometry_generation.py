"""Regression tests for deterministic STL geometry generation."""

from __future__ import annotations

import hashlib
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

from grid.stl_to_gmsh.generate_stl_geometry import generate_geometry
from grid.stl_to_gmsh.converter import convert


class StlGeometryGenerationTest(unittest.TestCase):
    def test_unit_sphere_matches_documented_icosphere_topology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = generate_geometry("sphere", Path(temporary) / "sphere.stl")
            self.assertEqual(result["surface"]["triangle_count"], 1280)
            self.assertEqual(result["surface"]["vertex_count"], 642)
            self.assertTrue(result["surface"]["watertight"])
            self.assertEqual(result["surface"]["connected_component_count"], 1)
            self.assertGreater(result["surface"]["signed_volume"], 0.0)
            self.assertEqual(result["patch_triangle_indices"]["wall"], list(range(1280)))
            self.assertEqual(result["units"], "m")
            self.assertEqual(result["stl_sha256"], hashlib.sha256(result["stl"].read_bytes()).hexdigest())
            self.assertTrue(result["manifest"].is_file())

    def test_prolate_spheroid_uses_configured_semi_axes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = generate_geometry(
                "prolate_spheroid",
                Path(temporary) / "prolate.stl",
                semi_axes=(3.0, 0.5, 0.5),
            )
            self.assertTrue(result["surface"]["watertight"])
            self.assertGreater(result["surface"]["bounds"]["x"][1], 2.85)
            self.assertLessEqual(result["surface"]["bounds"]["x"][1], 3.0)
            self.assertGreater(result["surface"]["bounds"]["y"][1], 0.475)
            self.assertLessEqual(result["surface"]["bounds"]["y"][1], 0.5)
            self.assertGreater(result["surface"]["bounds"]["z"][1], 0.475)
            self.assertLessEqual(result["surface"]["bounds"]["z"][1], 0.5)
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
            self.assertEqual(surface["connected_component_count"], 1)
            self.assertAlmostEqual(surface["signed_volume"], 32 * 0.5 * 0.5 * math.sin(2 * math.pi / 32) * 4.0 / 2, delta=0.02)
            self.assertEqual(len(labels["inflow"]), 32)
            self.assertEqual(len(labels["outflow"]), 32)
            self.assertEqual(len(labels["wall"]), 2 * 32 * 8)
            self.assertEqual(
                sum(map(len, labels.values())),
                surface["triangle_count"],
            )
            manifest = json.loads(result["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["intended_gmsh_physical_groups"]["surface"], ["wall", "inflow", "outflow"])
            self.assertIn("not solid walls", manifest["notes"])


    def test_cylinder_manifest_generates_distinct_gmsh_physical_surfaces(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stl = root / "cylinder.stl"
            generated = generate_geometry(
                "cylinder_internal_flow",
                stl,
                radius=0.5,
                length=4.0,
                circumferential_divisions=32,
                axial_divisions=8,
            )
            sidecar = root / "layout.json"
            sidecar.write_text(
                json.dumps({
                    "mesh_dimension": 3,
                    "volume_mesh": "cylinder_internal_flow",
                    "mesh_size_min": 0.1,
                    "mesh_size_max": 0.2,
                }),
                encoding="utf-8",
            )
            geo = root / "cylinder.geo"
            report = convert(stl, sidecar, geo)
            text = geo.read_text(encoding="utf-8")
            self.assertIn('Physical Surface("wall")', text)
            self.assertIn('Physical Surface("inflow")', text)
            self.assertIn('Physical Surface("outflow")', text)
            self.assertIn('Physical Volume("fluid_domain")', text)
            self.assertIn("Cylinder(1)", text)
            self.assertEqual(
                report["patch_triangle_counts"],
                {"wall": 512, "inflow": 32, "outflow": 32},
            )
            self.assertEqual(report["geometry"], "cylinder_internal_flow")
            self.assertTrue(os.path.samefile(generated["stl"], stl))

    def test_cylinder_manifest_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stl = root / "cylinder.stl"
            generate_geometry("cylinder_internal_flow", stl)
            stl.write_bytes(stl.read_bytes() + b"tampered")
            sidecar = root / "layout.json"
            sidecar.write_text(
                json.dumps({"mesh_dimension": 3, "volume_mesh": "cylinder_internal_flow"}),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "checksum"):
                convert(stl, sidecar, root / "cylinder.geo")

    def test_generated_stl_cannot_pollute_repository_source_tree(self) -> None:
        repository_root = Path(__file__).resolve().parents[4]
        protected_output = repository_root / "build" / "should_not_be_created.stl"
        with self.assertRaisesRegex(ValueError, "outside the OneFLOW source tree"):
            generate_geometry("sphere", protected_output)

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
