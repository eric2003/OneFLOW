"""Geometry and contract tests for the OECD/NEA-Vattenfall T-junction Case."""
from __future__ import annotations

import importlib.util
import json
import math
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
GENERATOR = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/python/generate_pipe_tee.py"
SPEC = importlib.util.spec_from_file_location("oneflow_pipe_tee_generator", GENERATOR)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PipeTeeGeometryTest(unittest.TestCase):
    def test_benchmark_dimensions_and_derived_bulk_velocities_are_preserved(self):
        p = MODULE.DEFAULTS
        self.assertAlmostEqual(p["main_diameter"], 0.140)
        self.assertAlmostEqual(p["branch_diameter"], 0.100)
        self.assertAlmostEqual(p["main_inlet_length"], 1.070)
        self.assertAlmostEqual(p["main_outlet_length"], 1.070)
        self.assertAlmostEqual(p["branch_inlet_length"], 0.470)
        main_area = math.pi * p["main_diameter"] ** 2 / 4.0
        branch_area = math.pi * p["branch_diameter"] ** 2 / 4.0
        self.assertAlmostEqual(0.009 / main_area, 0.5848, places=4)
        self.assertAlmostEqual(0.006 / branch_area, 0.7639, places=4)

    def test_geo_definition_has_boolean_union_and_three_named_surface_regions(self):
        text = MODULE._geo_text(MODULE.DEFAULTS)
        self.assertIn("BooleanUnion", text)
        self.assertIn('Physical Surface("wall"', text)
        self.assertIn('Physical Surface("inflow"', text)
        self.assertIn('Physical Surface("outflow"', text)
        self.assertIn("Surface In BoundingBox", text)

    def test_exporter_creates_closed_surface_and_checksum_bound_regions(self):
        from grid.stl_import.reader import read_stl
        from grid.stl_import.validation import analyze_surface

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            nodes = {
                1:(-1.0,-1.0,-1.0), 2:(1.0,-1.0,-1.0), 3:(1.0,1.0,-1.0), 4:(-1.0,1.0,-1.0),
                5:(-1.0,-1.0,1.0), 6:(1.0,-1.0,1.0), 7:(1.0,1.0,1.0), 8:(-1.0,1.0,1.0),
            }
            # A closed cube fixture, with one physical region per face group.
            groups = {
                "wall":[(1,2,6),(1,6,5),(2,3,7),(2,7,6),(3,4,8),(3,8,7),(4,1,5),(4,5,8)],
                "inflow":[(1,4,3),(1,3,2)],
                "outflow":[(5,6,7),(5,7,8)],
            }
            msh = root / "fixture.msh"
            names = [(2,1,"wall"),(2,2,"inflow"),(2,3,"outflow")]
            with msh.open("w",encoding="utf-8") as f:
                f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n3\n")
                for dim,tag,name in names: f.write(f'{dim} {tag} "{name}"\n')
                f.write("$EndPhysicalNames\n$Nodes\n8\n")
                for i,p in nodes.items(): f.write(f"{i} {p[0]} {p[1]} {p[2]}\n")
                f.write("$EndNodes\n$Elements\n12\n")
                index=1
                for tag,name in enumerate(("wall","inflow","outflow"),1):
                    for face in groups[name]:
                        f.write(f"{index} 2 2 {tag} {tag} {face[0]} {face[1]} {face[2]}\n")
                        index+=1
                f.write("$EndElements\n")
            output = root / "cube.stl"
            result = MODULE.export_surface(msh, output)
            self.assertTrue(output.is_file())
            self.assertEqual(set(result["surface"]["region_triangle_counts"]), {"wall","inflow","outflow"})
            manifest = json.loads((root / "region_manifest.json").read_text(encoding="utf-8"))
            for name, record in manifest["regions"].items():
                self.assertTrue((root / record["file"]).is_file())
                self.assertEqual(record["triangle_count"], len(read_stl(root / record["file"])))
            metrics = analyze_surface(read_stl(output))
            self.assertEqual(metrics["status"], "valid")
            self.assertTrue(metrics["closed"])
            self.assertTrue(metrics["edge_manifold"])
            self.assertEqual(metrics["connected_component_count"], 1)

    @unittest.skipUnless(shutil.which("gmsh"), "Gmsh is required for the geometry integration test")
    def test_gmsh_surface_mesh_contains_the_benchmark_physical_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            stl = Path(temporary) / "pipe_tee.stl"
            result = MODULE.generate(stl, shutil.which("gmsh"), mesh_size=0.045)
            self.assertTrue(stl.is_file())
            manifest = json.loads(result["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["asset_id"], "geometry.pipe_tee.oecd_nea_vattenfall")
            self.assertEqual(set(manifest["intended_gmsh_physical_groups"]["surface"]), {"wall","inflow","outflow"})


if __name__ == "__main__":
    unittest.main()
