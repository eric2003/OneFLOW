"""Generic plane cutting tests using the canonical cylinder STL as blind input."""
import contextlib
import io
import json
import math
import os
import shutil
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
        self.case_directory = (
            root / "grid" / "stl_to_gmsh" / "3d" / "cylinder_internal_flow"
            / "oblique_cut" / "python"
        )
        self.case_config = json.loads(
            (self.case_directory.parent / "script" / "cutting.json").read_text(encoding="utf-8")
        )

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


    def test_oblique_cylinder_cut_creates_elliptical_inlet_and_circular_outlet(self):
        config = self.case_config
        planes = config["cutting_planes"]
        inlet, outlet = planes["inlet"], planes["outlet"]
        result = cut_closed_surface_between_planes(
            self.source, tuple(inlet["point"]), tuple(inlet["normal"]),
            tuple(outlet["point"]), tuple(outlet["normal"]),
            tolerance=config["tolerance"],
            inlet_keep_positive=inlet.get("keep_positive", True),
            outlet_keep_positive=outlet.get("keep_positive", False),
        )
        angle = math.radians(inlet["inclination_degrees"])
        major_axis = (-math.sin(angle), math.cos(angle), 0.0)
        minor_axis = (0.0, 0.0, 1.0)
        center = tuple(inlet["point"])

        def semi_extent(loop, axis):
            values = [
                sum((point[i] - center[i]) * axis[i] for i in range(3))
                for point in loop
            ]
            return 0.5 * (max(values) - min(values))

        measured_ratio = (
            semi_extent(result["inlet_loop"], major_axis)
            / semi_extent(result["inlet_loop"], minor_axis)
        )
        expected_ratio = config["analytic_expectation_for_verification_only"][
            "inlet_major_to_minor_ratio"
        ]
        self.assertGreater(measured_ratio, 1.10)
        self.assertAlmostEqual(measured_ratio, expected_ratio, delta=0.04)

        outlet_y = [point[1] for point in result["outlet_loop"]]
        outlet_z = [point[2] for point in result["outlet_loop"]]
        outlet_ratio = (
            (max(outlet_y) - min(outlet_y))
            / (max(outlet_z) - min(outlet_z))
        )
        self.assertAlmostEqual(outlet_ratio, 1.0, delta=0.04)
        topology = validate_closed_genus_zero_surface(
            result["closed_triangles"], config["tolerance"]
        )
        self.assertEqual(topology["euler_characteristic"], 2)
        self.assertGreater(topology["signed_volume"], 0.0)
        self.assertGreater(result["inlet_cap_triangle_count"], 0)
        self.assertGreater(result["outlet_cap_triangle_count"], 0)

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
            self.assertIn("Coherence Mesh;", geo)
            self.assertIn("CreateTopology;", geo)
            self.assertLess(geo.index("CreateTopology;"), geo.index("CreateGeometry;"))
            self.assertIn('Physical Surface("wall")', geo)
            self.assertIn('Physical Surface("inflow")', geo)
            self.assertIn('Physical Surface("outflow")', geo)
            self.assertIn('Physical Volume("fluid_domain")', geo)
            self.assertEqual(
                set(geometry_report["region_triangle_counts"]),
                {"wall", "inflow", "outflow"},
            )


    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 and install Gmsh to run the volume-mesh integration test",
    )
    def test_gmsh_generates_volume_mesh_with_disjoint_named_boundaries(self):
        from grid.stl_to_gmsh.generate_volume_mesh import generate

        gmsh = shutil.which("gmsh") or shutil.which("gmsh.exe")
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

            mesh_path = output / "results" / "pipe.msh"
            report = generate(
                output / "closed_surface.stl",
                mesh_path,
                gmsh,
                mesh_size_min=0.20,
                mesh_size_max=0.30,
            )
            self.assertTrue(report["validated"], report)
            mesh = report["mesh_summary"]
            self.assertGreater(mesh["tetrahedron_elements"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            self.assertEqual(
                mesh["tetrahedral_topology"]["surface_triangles"],
                mesh["triangle_elements"],
            )

            groups = {
                (item["dimension"], item["name"]): item["element_count"]
                for item in mesh["physical_groups"]
            }
            expected = {
                (2, "wall"), (2, "inflow"), (2, "outflow"),
                (3, "fluid_domain"),
            }
            self.assertEqual(set(groups), expected)
            for group in expected:
                self.assertGreater(groups[group], 0, group)
            self.assertEqual(
                sum(groups[(2, name)] for name in ("wall", "inflow", "outflow")),
                mesh["triangle_elements"],
                "Every exterior triangle must belong to exactly one named boundary group",
            )


    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 and install Gmsh to run the curved-surface integration test",
    )
    def test_gmsh_meshes_blind_non_cylindrical_surface_with_two_cut_loops(self):
        from grid.stl_to_gmsh.generate_volume_mesh import generate

        root = Path(__file__).resolve().parents[4]
        source = (
            root / "grid" / "stl_to_gmsh" / "3d" / "prolate_spheroid"
            / "automatic_tetrahedra" / "python" / "input" / "model.stl"
        )
        gmsh = shutil.which("gmsh") or shutil.which("gmsh.exe")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "cut"
            arguments = [
                "--input", str(source),
                "--output-directory", str(output),
                "--inlet-point", "-2.5", "0", "0",
                "--inlet-normal", "1", "0", "0",
                "--outlet-point", "2.5", "0", "0",
                "--outlet-normal", "1", "0", "0",
            ]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cut_surface_main(arguments), 0)

            manifest = __import__("json").loads(
                (output / "region_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(set(manifest["regions"]), {"wall", "inflow", "outflow"})
            self.assertGreater(manifest["regions"]["wall"]["triangle_count"], 0)

            mesh_path = output / "results" / "curved_shell.msh"
            report = generate(
                output / "closed_surface.stl",
                mesh_path,
                gmsh,
                mesh_size_min=0.10,
                mesh_size_max=0.15,
            )
            self.assertTrue(report["validated"], report)
            mesh = report["mesh_summary"]
            self.assertGreater(mesh["tetrahedron_elements"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            groups = {
                (item["dimension"], item["name"]): item["element_count"]
                for item in mesh["physical_groups"]
            }
            expected = {
                (2, "wall"), (2, "inflow"), (2, "outflow"),
                (3, "fluid_domain"),
            }
            self.assertEqual(set(groups), expected)
            for group in expected:
                self.assertGreater(groups[group], 0, group)
            self.assertEqual(
                sum(groups[(2, name)] for name in ("wall", "inflow", "outflow")),
                mesh["triangle_elements"],
                "Every curved-shell exterior triangle must belong to one named boundary group",
            )


    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 and install Gmsh for the oblique-cut volume-mesh test",
    )
    def test_gmsh_meshes_oblique_cylinder_with_elliptical_inlet_group(self):
        from grid.stl_to_gmsh.generate_volume_mesh import generate

        config = self.case_config
        planes = config["cutting_planes"]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "oblique-cylinder"
            arguments = [
                "--input", str(self.case_directory / "input" / "model.stl"),
                "--output-directory", str(output),
                "--inlet-point", *(str(value) for value in planes["inlet"]["point"]),
                "--inlet-normal", *(str(value) for value in planes["inlet"]["normal"]),
                "--outlet-point", *(str(value) for value in planes["outlet"]["point"]),
                "--outlet-normal", *(str(value) for value in planes["outlet"]["normal"]),
                "--tolerance", str(config["tolerance"]),
            ]
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cut_surface_main(arguments), 0)
            report = generate(
                output / "closed_surface.stl",
                output / "results" / "cylinder_oblique_cut.msh",
                shutil.which("gmsh") or shutil.which("gmsh.exe"),
                mesh_size_min=config["volume_mesh"]["mesh_size_min"],
                mesh_size_max=config["volume_mesh"]["mesh_size_max"],
            )
            self.assertTrue(report["validated"], report)
            mesh = report["mesh_summary"]
            self.assertGreater(mesh["tetrahedron_elements"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            groups = {
                (item["dimension"], item["name"]): item["element_count"]
                for item in mesh["physical_groups"]
            }
            self.assertEqual(
                set(groups),
                {(2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")},
            )
            self.assertTrue(all(count > 0 for count in groups.values()))
            self.assertEqual(
                sum(groups[(2, name)] for name in ("wall", "inflow", "outflow")),
                mesh["triangle_elements"],
            )

    def test_empty_input_is_rejected(self):
        with self.assertRaises(ValueError):
            cut_closed_surface_between_planes(
                [], (0.,0.,0.), (1.,0.,0.), (1.,0.,0.), (1.,0.,0.)
            )


if __name__ == "__main__":
    unittest.main()
