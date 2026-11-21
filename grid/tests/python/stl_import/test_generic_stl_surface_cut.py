"""Generic plane cutting tests using the canonical cylinder STL as blind input."""
import unittest
from pathlib import Path

from grid.stl_import.reader import read_stl
from grid.stl_import.surface_cut import cut_closed_surface_between_planes
from grid.stl_to_gmsh.converter import validate_closed_genus_zero_surface


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

    def test_empty_input_is_rejected(self):
        with self.assertRaises(ValueError):
            cut_closed_surface_between_planes(
                [], (0.,0.,0.), (1.,0.,0.), (1.,0.,0.), (1.,0.,0.)
            )


if __name__ == "__main__":
    unittest.main()
