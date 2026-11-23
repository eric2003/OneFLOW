"""Tests for analytic cylinder volume-mesh generation contracts."""
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from grid.stl_to_gmsh.generate_cylinder_volume_mesh import (
    build_cylinder_volume_geo, _require_external_output, generate,
)


class CylinderVolumeMeshTest(unittest.TestCase):
    def test_geo_defines_cut_volume_and_explicit_physical_groups(self):
        geo = build_cylinder_volume_geo(1.0, -1.5, 1.5, 0.2)
        self.assertIn('Cylinder(1) = {-1.5, 0, 0, 3, 0, 0, 1};', geo)
        self.assertIn('Physical Surface("wall")', geo)
        self.assertIn('Physical Surface("inflow")', geo)
        self.assertIn('Physical Surface("outflow")', geo)
        self.assertIn('Physical Volume("fluid_domain")', geo)
        self.assertIn('Mesh.MeshSizeMax = 0.20000000000000001;', geo)

    def test_invalid_parameters_are_rejected(self):
        for values in ((0.0, -1.5, 1.5, 0.2), (1.0, 1.5, -1.5, 0.2),
                       (1.0, -1.5, 1.5, 0.0)):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    build_cylinder_volume_geo(*values)

    def test_output_inside_repository_is_rejected(self):
        repository_root = Path(__file__).resolve().parents[4]
        with self.assertRaises(ValueError):
            _require_external_output(repository_root / "grid" / "generated")


    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 to run the manual Gmsh integration test",
    )
    def test_gmsh_builds_volume_with_required_physical_groups(self):
        gmsh = shutil.which("gmsh") or shutil.which("gmsh.exe")
        with tempfile.TemporaryDirectory() as temporary:
            report = generate(
                Path(temporary) / "cylinder-cut",
                gmsh,
                cut_start=-1.5,
                cut_end=1.5,
                mesh_size=0.5,
                circumferential_divisions=24,
                axial_divisions=4,
            )
            self.assertEqual(report["stages"]["stage_a_boundary_loop_inspection"], "PASS")
            self.assertEqual(report["stages"]["gmsh_analytic_volume_mesh"], "PASS")
            self.assertGreater(report["mesh_summary"]["tetrahedron_elements"], 0)
            groups = {
                (item["dimension"], item["name"])
                for item in report["mesh_summary"]["physical_groups"]
            }
            self.assertEqual(groups, {
                (2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")
            })



if __name__ == "__main__":
    unittest.main()
