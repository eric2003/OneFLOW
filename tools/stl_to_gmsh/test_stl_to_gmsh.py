from __future__ import annotations
import json, struct, tempfile, unittest
from pathlib import Path
from .stl_to_gmsh import Stl2DError, build_geo, convert, extract_boundary, read_stl

DATA = Path(__file__).parent / "data"

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

if __name__ == "__main__": unittest.main()
