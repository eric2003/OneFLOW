"""Generic plane cutting tests using the canonical cylinder STL as blind input."""
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from grid.stl_import.reader import read_stl
from grid.stl_import.surface_cut import cut_closed_surface_between_planes
from grid.stl_import.cut_surface_cli import main as cut_surface_main
from grid.stl_to_gmsh.converter import (
    _region_aware_closed_surface_volume_geo,
    validate_closed_genus_zero_surface,
)


class GenericStlSurfaceCutTest(unittest.TestCase):
    def setUp(self):
        root = Path(__file__).resolve().parents[4]
        self.stl = (
            root / "grid" / "stl_to_gmsh" / "3d" / "cylinder_internal_flow"
            / "automatic_tetrahedra" / "python" / "input" / "model.stl"
        )
        self.source = read_stl(self.stl)

    def test_cylinder_stl_is_cut_without_reading_analytic_manifest(self):
        result = cut_closed_surface_between_planes(
            self.source, (-1.5,0.0,0.0), (1.0,0.0,0.0),
            (1.5,0.0,0.0), (1.0,0.0,0.0), tolerance=1.0e-8,
        )
        self.assertEqual(result["source_triangle_count"], len(self.source))
        self.assertEqual(result["inlet_loop_vertex_count"], 48)
        self.assertEqual(result["outlet_loop_vertex_count"], 48)
        self.assertEqual(result["inlet_cap_triangle_count"], 46)
        self.assertEqual(result["outlet_cap_triangle_count"], 46)
        topology = validate_closed_genus_zero_surface(result["closed_triangles"], 1.0e-8)
        self.assertEqual(topology["euler_characteristic"], 2)
        self.assertGreater(topology["signed_volume"], 0.0)

    def test_nonparallel_cut_planes_are_supported(self):
        result = cut_closed_surface_between_planes(
            self.source, (-1.5,0.0,0.0), (1.0,0.05,0.0),
            (1.5,0.0,0.0), (1.0,-0.05,0.0), tolerance=1.0e-8,
        )
        self.assertGreater(result["inlet_loop_vertex_count"], 3)
        self.assertGreater(result["outlet_loop_vertex_count"], 3)
        validate_closed_genus_zero_surface(result["closed_triangles"], 1.0e-8)

    def test_cli_writes_named_surface_regions_and_report(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "cut"
            arguments = [
                "--input", str(self.stl),
                "--output-directory", str(output),
                "--inlet-point", "-1.5", "0", "0",
                "--inlet-normal", "1", "0", "0",
                "--outlet-point", "1.5", "0", "0",
                "--outlet-normal", "1", "0", "0",
            ]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cut_surface_main(arguments), 0)
            for name in ("wall.stl", "inlet.stl", "outlet.stl", "closed_surface.stl", "cut_report.json"):
                self.assertTrue((output / name).is_file(), name)
            json = __import__("json")
            report = json.loads((output / "cut_report.json").read_text(encoding="utf-8"))
            self.assertFalse(report["analytic_geometry_metadata_read"])
            self.assertEqual(report["operation"], "generic_stl_surface_cut")
            region_manifest = json.loads(
                (output / "region_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                set(region_manifest["regions"]), {"wall", "inflow", "outflow"}
            )
            self.assertEqual(
                sum(item["triangle_count"] for item in region_manifest["regions"].values()),
                report["counts"]["side_triangle_count"]
                + report["counts"]["inlet_cap_triangle_count"]
                + report["counts"]["outlet_cap_triangle_count"],
            )
            closed_triangles = read_stl(output / "closed_surface.stl")
            geo, geometry_report = _region_aware_closed_surface_volume_geo(
                output / "closed_surface.stl",
                closed_triangles,
                {
                    "region_manifest": str(output / "region_manifest.json"),
                    "mesh_size_min": 0.1,
                    "mesh_size_max": 0.2,
                },
                1.0e-8,
            )
            self.assertIn('Physical Surface("wall")', geo)
            self.assertIn('Physical Surface("inflow")', geo)
            self.assertIn('Physical Surface("outflow")', geo)
            self.assertIn('Physical Volume("fluid_domain")', geo)
            self.assertEqual(
                set(geometry_report["region_triangle_counts"]),
                {"wall", "inflow", "outflow"},
            )

    def test_empty_input_is_rejected(self):
        with self.assertRaises(ValueError):
            cut_closed_surface_between_planes(
                [], (0.,0.,0.), (1.,0.,0.), (1.,0.,0.), (1.,0.,0.)
            )


if __name__ == "__main__":
    unittest.main()
