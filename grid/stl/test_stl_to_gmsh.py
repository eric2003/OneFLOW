from __future__ import annotations
import json, shutil, struct, subprocess, tempfile, unittest
from pathlib import Path
from .stl_to_gmsh import Stl2DError, _weld_vertices, build_geo, convert, extract_boundary, inspect_msh4, read_stl

DATA = Path(__file__).parent / "fixtures"

class Stl2DGmshTest(unittest.TestCase):
    def test_ascii_square_extracts_four_boundary_edges(self):
        tris = read_stl(DATA / "square.stl")
        self.assertEqual(len(tris), 2)
        self.assertEqual(len(extract_boundary(tris)), 4)

    def test_binary_stl_is_supported(self):
        tris = [((0.,0.,0.),(1.,0.,0.),(1.,1.,0.)),((0.,0.,0.),(1.,1.,0.),(0.,1.,0.))]
        data = bytearray(80) + bytearray(struct.pack("<I",2))
        for tri in tris: data.extend(struct.pack("<12fH",0,0,1,*tri[0],*tri[1],*tri[2],0))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "square.stl"; path.write_bytes(data)
            self.assertEqual(len(extract_boundary(read_stl(path))),4)

    def test_geo_has_quad_and_four_physical_boundary_groups(self):
        with tempfile.TemporaryDirectory() as temp:
            geo = Path(temp) / "square.geo"
            report = convert(DATA/"square.stl",DATA/"square.boundary.json",geo)
            content = geo.read_text()
            self.assertEqual(report["boundary_edges"],4)
            for name in ("Inlet","Outlet","Symmetry","Farfield"):
                self.assertIn(f'Physical Curve("{name}")',content)
            self.assertIn("Transfinite Surface {1};",content)
            self.assertIn("Recombine Surface {1};",content)
            self.assertIn('Physical Surface("Fluid")',content)

    def test_gmsh_cli_generates_one_quad_and_four_boundary_lines(self):
        gmsh = shutil.which("gmsh")
        if not gmsh:
            self.skipTest("Gmsh executable is not installed")
        with tempfile.TemporaryDirectory() as temp:
            geo = Path(temp) / "square.geo"
            msh = Path(temp) / "square.msh"
            convert(DATA/"square.stl", DATA/"square.boundary.json", geo)
            subprocess.run([gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)], check=True)
            summary = inspect_msh4(msh)
            self.assertEqual(summary["format"], "4.1")
            self.assertEqual(summary["nodes"], 4)
            self.assertEqual(summary["quadrilateral_elements"], 1)
            self.assertEqual(summary["line_elements"], 4)
            names = {group["name"] for group in summary["physical_groups"]}
            self.assertEqual(names, {"Inlet", "Outlet", "Symmetry", "Farfield", "Fluid"})

    def test_two_quad_layout_splits_the_horizontal_edges_into_two_cells(self):
        config = json.loads((DATA / "square.boundary.json").read_text(encoding="utf-8"))
        config.pop("one_quad_cell", None)
        config["transfinite_cells"] = {"x": 2, "y": 1}
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "two_quads.json"
            geo = Path(temp) / "two_quads.geo"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            report = convert(DATA / "square.stl", sidecar, geo)
            content = geo.read_text(encoding="utf-8")
            self.assertEqual(report["structured_cells"], {"x": 2, "y": 1, "total": 2})
            self.assertIn("Transfinite Curve {1, 3} = 3;", content)
            self.assertIn("Transfinite Curve {2, 4} = 2;", content)
            self.assertIn("Transfinite Surface {1};", content)
            self.assertIn("Recombine Surface {1};", content)

    def test_transfinite_cell_counts_must_be_positive_integers(self):
        config = json.loads((DATA / "square.boundary.json").read_text(encoding="utf-8"))
        config.pop("one_quad_cell", None)
        config["transfinite_cells"] = {"x": 0, "y": 1}
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "invalid.json"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(Stl2DError, "positive integers"):
                convert(DATA / "square.stl", sidecar, Path(temp) / "invalid.geo")

    def test_gmsh_cli_generates_two_adjacent_square_cells_and_boundary_groups(self):
        gmsh = shutil.which("gmsh")
        if not gmsh:
            self.skipTest("Gmsh executable is not installed")
        case = Path(__file__).parent / "stl_to_gmsh"
        with tempfile.TemporaryDirectory() as temp:
            geo = Path(temp) / "two_quads.geo"
            msh = Path(temp) / "two_quads.msh"
            convert(
                case / "input" / "two_square_cells.stl",
                case / "script" / "two_quads.json",
                geo,
            )
            subprocess.run([gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)], check=True)
            summary = inspect_msh4(msh)
            self.assertEqual(summary["format"], "4.1")
            self.assertEqual(summary["nodes"], 6)
            self.assertEqual(summary["line_elements"], 6)
            self.assertEqual(summary["quadrilateral_elements"], 2)
            names = {group["name"] for group in summary["physical_groups"]}
            self.assertEqual(names, {"inflow", "outflow", "symmetry", "farfield", "fluid_domain"})

    def test_stl_to_gmsh_case_uses_two_adjacent_unit_square_cells(self):
        case = Path(__file__).parent / "stl_to_gmsh"
        triangles = read_stl(case / "input" / "two_square_cells.stl")
        loop = extract_boundary(triangles)
        self.assertEqual(len(loop), 4)
        self.assertEqual((min(p[0] for p in loop), max(p[0] for p in loop)), (0.0, 2.0))
        self.assertEqual((min(p[1] for p in loop), max(p[1] for p in loop)), (0.0, 1.0))
        with tempfile.TemporaryDirectory() as temp:
            report = convert(
                case / "input" / "two_square_cells.stl",
                case / "script" / "two_quads.json",
                Path(temp) / "two_quads.geo",
            )
        self.assertEqual(report["structured_cells"], {"x": 2, "y": 1, "total": 2})

    def test_inspect_msh4_reads_element_blocks_and_boundary_names(self):
        content = """$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
1 1 "inflow"
1 2 "outflow"
1 3 "symmetry"
1 4 "farfield"
2 5 "fluid_domain"
$EndPhysicalNames
$Nodes
4 4 1 4
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
$EndNodes
$Elements
5 5 1 5
1 1 1 1
1 1 2
1 2 1 1
2 2 3
1 3 1 1
3 3 4
1 4 1 1
4 4 1
2 1 3 1
5 1 2 3 4
$EndElements
"""
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "square.msh"
            path.write_text(content, encoding="utf-8")
            summary = inspect_msh4(path)
        self.assertEqual(summary["nodes"], 4)
        self.assertEqual(summary["line_elements"], 4)
        self.assertEqual(summary["quadrilateral_elements"], 1)
        self.assertEqual(
            {group["name"] for group in summary["physical_groups"]},
            {"inflow", "outflow", "symmetry", "farfield", "fluid_domain"},
        )

    def test_non_planar_stl_is_rejected(self):
        with self.assertRaisesRegex(Stl2DError,"constant-Z"):
            extract_boundary([((0.,0.,0.),(1.,0.,0.),(1.,1.,0.1))])

    def test_disconnected_loops_are_rejected(self):
        tris=[((0.,0.,0.),(1.,0.,0.),(0.,1.,0.)),((2.,0.,0.),(3.,0.,0.),(2.,1.,0.))]
        with self.assertRaises(Stl2DError): extract_boundary(tris)

    def test_sidecar_must_assign_every_edge(self):
        cfg=json.loads((DATA/"square.boundary.json").read_text()); cfg["boundaries"].pop()
        with tempfile.TemporaryDirectory() as temp:
            sidecar=Path(temp)/"partial.json"; sidecar.write_text(json.dumps(cfg))
            with self.assertRaisesRegex(Stl2DError,"every outer edge"):
                convert(DATA/"square.stl",sidecar,Path(temp)/"out.geo")

    def test_one_quad_requires_four_vertices(self):
        loop=[(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0.5,1.5,0.),(0.,1.,0.)]
        with self.assertRaisesRegex(Stl2DError,"exactly four"):
            build_geo(loop,{},{"one_quad_cell":True})


    def test_stl_body_generates_outer_fluid_domain_and_five_physical_groups(self):
        config = {
            "surface_name": "Fluid",
            "mesh_size": 0.25,
            "domain": {
                "xmin": -1.0, "xmax": 2.0, "ymin": -1.0, "ymax": 2.0,
                "boundaries": {
                    "left": "Inlet",
                    "right": "Outlet",
                    "bottom": "Symmetry",
                    "top": "Farfield"
                }
            },
            "body_boundary_name": "Wall"
        }
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "domain.json"
            geo = Path(temp) / "body_domain.geo"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            report = convert(DATA / "square.stl", sidecar, geo)
            content = geo.read_text(encoding="utf-8")
            self.assertEqual(report["body_boundary_edges"], 4)
            self.assertIn("Curve Loop(1) = {5, 6, 7, 8};", content)
            self.assertIn("Curve Loop(2) = {1, 2, 3, 4};", content)
            self.assertIn("Plane Surface(1) = {1, 2};", content)
            for name in ("Wall", "Inlet", "Outlet", "Symmetry", "Farfield", "Fluid"):
                self.assertIn(f'Physical Curve("{name}")' if name != "Fluid"
                              else f'Physical Surface("{name}")', content)

    def test_outer_domain_must_strictly_enclose_stl_body(self):
        config = {
            "domain": {
                "xmin": 0.0, "xmax": 2.0, "ymin": -1.0, "ymax": 2.0,
                "boundaries": {
                    "left": "Inlet", "right": "Outlet",
                    "bottom": "Symmetry", "top": "Farfield"
                }
            }
        }
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "domain.json"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(Stl2DError, "strictly enclose"):
                convert(DATA / "square.stl", sidecar, Path(temp) / "invalid.geo")

    def test_outer_domain_requires_all_four_distinct_boundary_names(self):
        config = {
            "domain": {
                "xmin": -1.0, "xmax": 2.0, "ymin": -1.0, "ymax": 2.0,
                "boundaries": {
                    "left": "Inlet", "right": "Outlet",
                    "bottom": "Farfield", "top": "Farfield"
                }
            }
        }
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "domain.json"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(Stl2DError, "four distinct"):
                convert(DATA / "square.stl", sidecar, Path(temp) / "invalid.geo")


    def test_gmsh_cli_generates_external_flow_domain_physical_groups(self):
        gmsh = shutil.which("gmsh")
        if not gmsh:
            self.skipTest("Gmsh executable is not installed")
        config = {
            "mesh_size": 0.5,
            "domain": {
                "xmin": -1.0, "xmax": 2.0, "ymin": -1.0, "ymax": 2.0,
                "boundaries": {
                    "left": "Inlet", "right": "Outlet",
                    "bottom": "Symmetry", "top": "Farfield"
                }
            },
            "body_boundary_name": "Wall"
        }
        with tempfile.TemporaryDirectory() as temp:
            sidecar = Path(temp) / "domain.json"
            geo = Path(temp) / "body_domain.geo"
            msh = Path(temp) / "body_domain.msh"
            sidecar.write_text(json.dumps(config), encoding="utf-8")
            convert(DATA / "square.stl", sidecar, geo)
            subprocess.run([gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)], check=True)
            summary = inspect_msh4(msh)
            names = {group["name"] for group in summary["physical_groups"]}
            self.assertEqual(names, {"Wall", "Inlet", "Outlet", "Symmetry", "Farfield", "Fluid"})
            self.assertGreater(summary["line_elements"], 0)
            self.assertGreater(summary["triangle_elements"], 0)

    def test_stl_inspection_case_writes_report_without_touching_source_tree(self):
        import sys
        script = Path(__file__).parent / "script" / "inspect.py"
        with tempfile.TemporaryDirectory() as temp:
            report = Path(temp) / "stl_report.json"
            subprocess.run(
                [sys.executable, str(script), "--input", str(DATA / "square.stl"),
                 "--output", str(report)],
                check=True,
            )
            result = json.loads(report.read_text(encoding="utf-8"))
            self.assertEqual(result["format"], "ascii")
            self.assertEqual(result["triangle_count"], 2)
            self.assertEqual(result["boundary_edge_count"], 4)


    def test_welding_matches_across_rounding_cell_boundary(self):
        raw = [
            (0.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.2, 0.2, 0.0),
            (0.201, 0.2, 0.0), (0.0, 0.2, 0.0), (0.001, 0.0, 0.0),
        ]
        points, triangles = _weld_vertices(raw, 0.1)
        self.assertEqual(len(points), 4)
        self.assertEqual(triangles[0][0], triangles[1][2])
        self.assertEqual(triangles[0][2], triangles[1][0])

    def test_welding_rejects_false_merge_within_rounded_cell(self):
        raw = [
            (0.49, 0.49, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
            (-0.49, -0.49, 0.0), (2.0, 1.0, 0.0), (0.0, 2.0, 1.0),
        ]
        points, triangles = _weld_vertices(raw, 1.0)
        self.assertNotEqual(triangles[0][0], triangles[1][0])
        self.assertEqual(len(points), 6)

    def test_welding_does_not_merge_at_exact_tolerance(self):
        points, triangles = _weld_vertices(
            [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 3.0, 0.0),
             (1.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0)],
            1.0,
        )
        self.assertNotEqual(triangles[0][0], triangles[1][0])

    def test_welding_does_not_take_transitive_closure(self):
        raw = [
            (0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 4.0, 0.0),
            (0.75, 0.0, 0.0), (4.0, 1.0, 0.0), (0.0, 4.0, 1.0),
            (1.5, 0.0, 0.0), (5.0, 0.0, 0.0), (0.0, 5.0, 0.0),
        ]
        points, triangles = _weld_vertices(raw, 1.0)
        self.assertEqual(triangles[0][0], triangles[1][0])
        self.assertNotEqual(triangles[0][0], triangles[2][0])

if __name__ == "__main__": unittest.main()
