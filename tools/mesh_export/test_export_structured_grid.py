import struct
import tempfile
import unittest
from pathlib import Path

from export_structured_grid import (
    build_topology,
    coordinate_bounds,
    read_structured_grid,
    write_svg,
)


class StructuredGridExportTest(unittest.TestCase):
    def test_minimal_rectangle_topology_counts(self):
        cells, faces, boundary_faces = build_topology(3, 3)

        self.assertEqual(len(cells), 4)
        self.assertEqual(len(faces), 12)
        self.assertEqual(len(boundary_faces), 8)
        self.assertEqual(sum(face["neighbour"] != -1 for face in faces), 4)
        self.assertTrue(all(len(cell["nodes"]) == 4 for cell in cells))

    def test_coordinate_bounds_include_all_three_axes(self):
        points = [
            [-1.0, 2.0, 0.0],
            [3.0, -4.0, 0.0],
            [0.5, 1.0, 0.0],
        ]

        self.assertEqual(
            coordinate_bounds(points),
            {"min": [-1.0, -4.0, 0.0], "max": [3.0, 2.0, 0.0]},
        )

    def test_empty_coordinate_bounds_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "empty mesh"):
            coordinate_bounds([])

    def test_svg_contains_node_markers_and_cell_ids(self):
        points = [
            [0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0],
            [0.0, 0.5, 0.0], [0.5, 0.5, 0.0], [1.0, 0.5, 0.0],
            [0.0, 1.0, 0.0], [0.5, 1.0, 0.0], [1.0, 1.0, 0.0],
        ]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "mesh.svg"
            write_svg(output, 3, 3, points)
            svg = output.read_text(encoding="utf-8")

        self.assertEqual(svg.count("<circle "), 9)
        self.assertEqual(svg.count("<text "), 4)
        self.assertIn('data-node-id="8"', svg)
        self.assertIn('data-cell-id="3"', svg)

    def test_reads_minimal_rectangle_binary_fixture(self):
        ni = nj = 3
        nk = 1
        x = [0.0, 0.5, 1.0] * 3
        y = [0.0] * 3 + [0.5] * 3 + [1.0] * 3
        z = [0.0] * 9
        payload = struct.pack("<4i", 1, ni, nj, nk)
        payload += struct.pack("<9d", *x)
        payload += struct.pack("<9d", *y)
        payload += struct.pack("<9d", *z)

        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "minimalRectangle.ofl"
            source.write_bytes(payload)
            actual_ni, actual_nj, points = read_structured_grid(source, "<", 8)

        self.assertEqual((actual_ni, actual_nj), (3, 3))
        self.assertEqual(len(points), 9)
        self.assertEqual(points[0], [0.0, 0.0, 0.0])
        self.assertEqual(points[-1], [1.0, 1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
