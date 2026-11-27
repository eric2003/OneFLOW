"""Geometry and contract tests for the OECD/NEA-Vattenfall T-junction Case."""
from __future__ import annotations

import importlib.util
import json
import math
import os
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
        self.assertAlmostEqual(p["main_inlet_length"], 1.220)
        self.assertAlmostEqual(p["main_outlet_length"], 1.220)
        self.assertAlmostEqual(p["branch_inlet_length"], 0.693)
        benchmark = json.loads((ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/automatic_tetrahedra/benchmark.json").read_text(encoding="utf-8"))
        self.assertEqual(p["main_diameter"], benchmark["geometry"]["main_pipe"]["inner_diameter"])
        self.assertEqual(p["branch_diameter"], benchmark["geometry"]["branch_pipe"]["inner_diameter"])
        self.assertEqual(p["main_inlet_length"], benchmark["geometry"]["main_pipe"]["inlet_length"])
        self.assertEqual(p["main_outlet_length"], benchmark["geometry"]["main_pipe"]["outlet_length"])
        self.assertEqual(p["branch_inlet_length"], benchmark["geometry"]["branch_pipe"]["inlet_length"])
        self.assertEqual(MODULE._BENCHMARK["id"], "oecd_nea_vattenfall_t_junction")
        self.assertAlmostEqual(p["mesh_size"], 0.030)
        main_area = math.pi * p["main_diameter"] ** 2 / 4.0
        branch_area = math.pi * p["branch_diameter"] ** 2 / 4.0
        self.assertAlmostEqual(0.009 / main_area, 0.5847, places=4)
        self.assertAlmostEqual(0.006 / branch_area, 0.7639, places=4)
        self.assertAlmostEqual(0.015 / main_area, 0.9744, places=4)

    def test_geo_definition_has_boolean_union_and_separate_inlet_surface_regions(self):
        text = MODULE._geo_text(MODULE.DEFAULTS)
        self.assertIn("Cylinder(1) = {-1.22, 0, 0, 2.44, 0, 0, 0.07};", text)
        self.assertIn("Cylinder(2) = {0, 0, 0, 0, 0, 0.693, 0.05};", text)
        self.assertIn("BooleanUnion", text)
        self.assertIn('Physical Surface("wall"', text)
        self.assertIn('Physical Surface("main_inflow"', text)
        self.assertIn('Physical Surface("branch_inflow"', text)
        self.assertIn('Physical Surface("outflow"', text)
        self.assertIn("Surface In BoundingBox", text)

    def test_surface_geometry_validator_checks_three_orthogonal_ports(self):
        nodes = {}
        groups = {"wall": [], "main_inflow": [], "branch_inflow": [], "outflow": []}
        params = dict(MODULE.DEFAULTS)
        next_id = 1

        def add_cap(name, axis, coordinate, radius):
            nonlocal next_id
            center = [0.0, 0.0, 0.0]
            center[axis] = coordinate
            center_id = next_id
            nodes[center_id] = tuple(center)
            next_id += 1
            radial_axes = [index for index in range(3) if index != axis]
            ring = []
            for index in range(64):
                angle = 2.0 * math.pi * index / 64
                point = [0.0, 0.0, 0.0]
                point[axis] = coordinate
                point[radial_axes[0]] = radius * math.cos(angle)
                point[radial_axes[1]] = radius * math.sin(angle)
                ring.append(next_id)
                nodes[next_id] = tuple(point)
                next_id += 1
            for index in range(64):
                groups[name].append((center_id, ring[index], ring[(index + 1) % 64]))
            return ring

        add_cap("main_inflow", 0, -params["main_inlet_length"], params["main_diameter"] / 2)
        add_cap("outflow", 0, params["main_outlet_length"], params["main_diameter"] / 2)
        branch_ring = add_cap("branch_inflow", 2, params["branch_inlet_length"], params["branch_diameter"] / 2)
        def add_wall_face(points):
            nonlocal next_id
            face = []
            for point in points:
                face.append(next_id)
                nodes[next_id] = point
                next_id += 1
            groups["wall"].append(tuple(face))

        add_wall_face(((-0.50, 0.070, 0.0), (-0.49, 0.070 * math.cos(0.1), 0.070 * math.sin(0.1)), (-0.48, 0.070 * math.cos(0.2), 0.070 * math.sin(0.2))))
        add_wall_face(((0.050, 0.0, 0.15), (0.050 * math.cos(0.2), 0.050 * math.sin(0.2), 0.15), (0.050 * math.cos(0.4), 0.050 * math.sin(0.4), 0.15)))

        report = MODULE._validate_pipe_tee_surface_geometry(nodes, groups, params)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["axis_angle_degrees"], 90.0)
        self.assertEqual(report["port_count"], 3)
        self.assertEqual(report["ports"]["branch_inflow"]["axis"], "+Z")

        # A branch cap moved onto a non-Z plane must be rejected.
        for node_id in branch_ring:
            x, y, z = nodes[node_id]
            nodes[node_id] = (params["branch_inlet_length"], y, 0.0)
        with self.assertRaisesRegex(ValueError, "branch_inflow.*end plane"):
            MODULE._validate_pipe_tee_surface_geometry(nodes, groups, params)

    def test_surface_geometry_validator_rejects_internal_overlap_wall(self):
        nodes = {}
        groups = {"wall": [], "main_inflow": [], "branch_inflow": [], "outflow": []}
        params = dict(MODULE.DEFAULTS)
        next_id = 1

        def add_cap(name, axis, coordinate, radius):
            nonlocal next_id
            center = [0.0, 0.0, 0.0]
            center[axis] = coordinate
            center_id = next_id
            nodes[next_id] = tuple(center)
            next_id += 1
            radial_axes = [index for index in range(3) if index != axis]
            ring = []
            for index in range(64):
                angle = 2.0 * math.pi * index / 64
                point = [0.0, 0.0, 0.0]
                point[axis] = coordinate
                point[radial_axes[0]] = radius * math.cos(angle)
                point[radial_axes[1]] = radius * math.sin(angle)
                ring.append(next_id)
                nodes[next_id] = tuple(point)
                next_id += 1
            for index in range(64):
                groups[name].append((center_id, ring[index], ring[(index + 1) % 64]))

        add_cap("main_inflow", 0, -params["main_inlet_length"], params["main_diameter"] / 2)
        add_cap("outflow", 0, params["main_outlet_length"], params["main_diameter"] / 2)
        add_cap("branch_inflow", 2, params["branch_inlet_length"], params["branch_diameter"] / 2)
        face = []
        # These points lie on the main-pipe cylinder but inside the branch bore.
        for point in ((0.0, 0.0, 0.070), (0.001, 0.070 * math.sin(0.01), 0.070 * math.cos(0.01)), (-0.001, 0.070 * math.sin(0.02), 0.070 * math.cos(0.02))):
            face.append(next_id)
            nodes[next_id] = point
            next_id += 1
        groups["wall"] = [tuple(face)]
        with self.assertRaisesRegex(ValueError, "inside the branch bore"):
            MODULE._validate_pipe_tee_surface_geometry(nodes, groups, params)

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
                "main_inflow":[(1,4,3)],
                "branch_inflow":[(1,3,2)],
                "outflow":[(5,6,7),(5,7,8)],
            }
            msh = root / "fixture.msh"
            names = [(2,1,"wall"),(2,2,"main_inflow"),(2,3,"branch_inflow"),(2,4,"outflow")]
            with msh.open("w",encoding="utf-8") as f:
                f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n4\n")
                for dim,tag,name in names: f.write(f'{dim} {tag} "{name}"\n')
                f.write("$EndPhysicalNames\n$Nodes\n8\n")
                for i,p in nodes.items(): f.write(f"{i} {p[0]} {p[1]} {p[2]}\n")
                f.write("$EndNodes\n$Elements\n12\n")
                index=1
                for tag,name in enumerate(("wall","main_inflow","branch_inflow","outflow"),1):
                    for face in groups[name]:
                        f.write(f"{index} 2 2 {tag} {tag} {face[0]} {face[1]} {face[2]}\n")
                        index+=1
                f.write("$EndElements\n")
            output = root / "cube.stl"
            result = MODULE.export_surface(msh, output)
            self.assertTrue(output.is_file())
            self.assertEqual(set(result["surface"]["region_triangle_counts"]), {"wall", "inflow", "outflow"})
            manifest = json.loads((root / "region_manifest.json").read_text(encoding="utf-8"))
            for name, record in manifest["regions"].items():
                self.assertTrue((root / record["file"]).is_file())
                self.assertEqual(record["triangle_count"], len(read_stl(root / record["file"])))
            self.assertEqual(set(manifest["regions"]), {"wall", "inflow", "outflow"})
            self.assertEqual(set(manifest["subregions"]), {"main_inflow", "branch_inflow"})
            self.assertEqual({record["parent_region"] for record in manifest["subregions"].values()}, {"inflow"})
            self.assertEqual(manifest["regions"]["inflow"]["triangle_count"], 2)
            self.assertEqual(manifest["subregions"]["main_inflow"]["triangle_count"], 1)
            self.assertEqual(manifest["subregions"]["branch_inflow"]["triangle_count"], 1)
            for name, record in manifest["subregions"].items():
                self.assertTrue((root / record["file"]).is_file())
                self.assertEqual(record["triangle_count"], len(read_stl(root / record["file"])))
            metrics = analyze_surface(read_stl(output))
            geometry_manifest = json.loads(output.with_suffix(".geometry.json").read_text(encoding="utf-8"))
            self.assertEqual(geometry_manifest["mesh_generation"]["surface_mesh_size_m"], MODULE.DEFAULTS["mesh_size"])
            self.assertEqual(geometry_manifest["generator"]["kind"], "python")
            self.assertEqual(metrics["status"], "valid")
            self.assertTrue(metrics["closed"])
            self.assertTrue(metrics["edge_manifold"])
            self.assertTrue(metrics["consistently_oriented"])
            self.assertEqual(metrics["connected_component_count"], 1)

    def test_cpp_geometry_driver_compiles(self):
        candidates = ("cl", "c++", "g++", "clang++") if os.name == "nt" else ("c++", "g++", "clang++")
        compiler = next((shutil.which(name) for name in candidates if shutil.which(name)), None)
        if not compiler:
            self.skipTest("A C++17 compiler is required.")
        source = ROOT / "grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/geometry/cpp/generate_pipe_tee.cpp"
        with tempfile.TemporaryDirectory() as temporary:
            is_msvc = Path(compiler).name.lower() in {"cl", "cl.exe"}
            executable = Path(temporary) / ("pipe_tee_generator.exe" if os.name == "nt" else "pipe_tee_generator")
            if is_msvc:
                command = [compiler, "/nologo", "/std:c++17", "/EHsc", str(source), f"/Fe:{executable}"]
            else:
                command = [compiler, "-std=c++17", str(source), "-o", str(executable)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertTrue(executable.is_file())

    @unittest.skipUnless(shutil.which("gmsh"), "Gmsh is required for the geometry integration test")
    def test_gmsh_surface_mesh_contains_the_benchmark_physical_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            stl = Path(temporary) / "pipe_tee.stl"
            result = MODULE.generate(stl, shutil.which("gmsh"), mesh_size=0.045)
            self.assertTrue(stl.is_file())
            manifest = json.loads(result["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["asset_id"], "geometry.pipe_tee.oecd_nea_vattenfall")
            self.assertEqual(manifest["geometry_validation"]["status"], "PASS")
            self.assertEqual(manifest["geometry_validation"]["axis_angle_degrees"], 90.0)
            self.assertEqual(manifest["geometry_validation"]["internal_overlap_wall_triangles"], 0)
            self.assertGreater(manifest["geometry_validation"]["wall_surface_triangle_counts"]["main_pipe"], 0)
            self.assertGreater(manifest["geometry_validation"]["wall_surface_triangle_counts"]["branch_pipe"], 0)
            self.assertEqual(set(manifest["geometry_validation"]["ports"]), {"main_inflow", "branch_inflow", "outflow"})
            self.assertEqual(set(manifest["intended_gmsh_physical_groups"]["surface"]), {"wall", "main_inflow", "branch_inflow", "outflow"})
            self.assertEqual(set(manifest["intended_gmsh_physical_groups"]["direct_occ_volume_mesh_surface_groups"]), {"wall", "main_inflow", "branch_inflow", "outflow"})
            self.assertTrue(all(count > 0 for count in manifest["surface"]["inlet_subregion_triangle_counts"].values()))
            from grid.stl_to_gmsh.converter import inspect_msh4
            geo = stl.with_suffix(".geo")
            volume_path = Path(temporary) / "pipe_tee.msh"
            volume_geo = Path(temporary) / "pipe_tee.volume.geo"
            volume_geo.write_text(geo.read_text(encoding="utf-8") + "\nMesh.MshFileVersion = 4.1;\n",
                                  encoding="utf-8")
            subprocess.run([shutil.which("gmsh"), str(volume_geo), "-3", "-format", "msh4", "-o", str(volume_path)],
                           check=True, capture_output=True, text=True)
            benchmark_geometry = MODULE._BENCHMARK["geometry"]
            volume = inspect_msh4(volume_path, geometry_contract={
                "main_diameter": benchmark_geometry["main_pipe"]["inner_diameter"],
                "branch_diameter": benchmark_geometry["branch_pipe"]["inner_diameter"],
                "main_inlet_length": benchmark_geometry["main_pipe"]["inlet_length"],
                "main_outlet_length": benchmark_geometry["main_pipe"]["outlet_length"],
                "branch_inlet_length": benchmark_geometry["branch_pipe"]["inlet_length"],
            })
            groups = {(item["dimension"], item["name"]) for item in volume["physical_groups"]}
            self.assertEqual(groups, {(2,"wall"),(2,"main_inflow"),(2,"branch_inflow"),(2,"outflow"),(3,"fluid_domain")})
            self.assertGreater(volume["tetrahedron_elements"], 0)
            self.assertEqual(volume["geometry_validation"]["status"], "PASS")
            self.assertGreater(volume["geometry_validation"]["checked_tetrahedra"], 0)
            self.assertIn('Physical Volume("fluid_domain"', geo.read_text(encoding="utf-8"))

            # Exercise both UI operating modes with the real Gmsh parser and
            # verify that the two main ports receive one shared Physical Group.
            import runpy
            runner_path = (
                Path(__file__).resolve().parents[4] / "grid" / "stl_to_gmsh" / "3d"
                / "pipe_tee" / "oecd_nea_vattenfall" / "automatic_tetrahedra"
                / "run_pipe_tee_case.py"
            )
            configure_mode = runpy.run_path(str(runner_path))["configure_operating_mode"]
            expected_mode_groups = {
                "main_inlets": {
                    (2, "wall"), (2, "main_left_inflow"), (2, "main_right_inflow"),
                    (2, "branch_outflow"), (3, "fluid_domain")
                },
                "branch_inlet": {
                    (2, "wall"), (2, "main_left_outflow"), (2, "main_right_outflow"),
                    (2, "branch_inflow"), (3, "fluid_domain")
                },
            }
            for mode, expected in expected_mode_groups.items():
                mode_geo = Path(temporary) / f"pipe_tee_{mode}.geo"
                mode_msh = Path(temporary) / f"pipe_tee_{mode}.msh"
                mode_geo.write_text(
                    geo.read_text(encoding="utf-8"), encoding="utf-8"
                )
                configure_mode(mode_geo, mode)
                mode_volume_geo = mode_geo.with_name(f"pipe_tee_{mode}.volume.geo")
                mode_volume_geo.write_text(
                    mode_geo.read_text(encoding="utf-8") + "\nMesh.MshFileVersion = 4.1;\n",
                    encoding="utf-8",
                )
                subprocess.run(
                    [shutil.which("gmsh"), str(mode_volume_geo), "-3", "-format", "msh4", "-o", str(mode_msh)],
                    check=True, capture_output=True, text=True,
                )
                mode_volume = inspect_msh4(mode_msh, geometry_contract={
                    "main_diameter": benchmark_geometry["main_pipe"]["inner_diameter"],
                    "branch_diameter": benchmark_geometry["branch_pipe"]["inner_diameter"],
                    "main_inlet_length": benchmark_geometry["main_pipe"]["inlet_length"],
                    "main_outlet_length": benchmark_geometry["main_pipe"]["outlet_length"],
                    "branch_inlet_length": benchmark_geometry["branch_pipe"]["inlet_length"],
                })
                mode_groups = {
                    (item["dimension"], item["name"])
                    for item in mode_volume["physical_groups"]
                }
                self.assertEqual(mode_groups, expected)
                self.assertGreater(mode_volume["tetrahedron_elements"], 0)
                self.assertEqual(mode_volume["geometry_validation"]["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
