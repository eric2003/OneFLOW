"""Tests for explicit layout topology validation."""

from __future__ import annotations

import unittest

from tools.mesh_layout.layout_topology import (
    LayoutBlock,
    LayoutFace,
    LayoutLine,
    LayoutPoint,
    LayoutTopology,
)


class LayoutTopologyTest(unittest.TestCase):
    def make_topology(self):
        return LayoutTopology(
            points=(
                LayoutPoint(1, 0.0, 0.0, 0.0),
                LayoutPoint(2, 1.0, 0.0, 0.0),
                LayoutPoint(3, 0.0, 1.0, 0.0),
            ),
            lines=(
                LayoutLine(1, 1, 2, 3, ("auto", 0.1), 2),
                LayoutLine(2, 2, 3, 3, ("auto", 0.1), 2),
                LayoutLine(3, 3, 1, 3, ("auto", 0.1), 2),
            ),
            faces=(LayoutFace(1, ((1, 1), (2, 2), (3, 3))),),
            blocks=(LayoutBlock(1, ((1, 1),)),),
        )

    def test_valid_point_line_face_block_references(self):
        self.make_topology().validate()

    def test_line_must_reference_existing_points(self):
        topology = self.make_topology()
        invalid_lines = (
            LayoutLine(1, 1, 9, 3, (), 2),
            topology.lines[1],
            topology.lines[2],
        )
        invalid = LayoutTopology(
            topology.points, invalid_lines, topology.faces, topology.blocks
        )
        with self.assertRaisesRegex(ValueError, "unknown point"):
            invalid.validate()

    def test_face_must_reference_existing_lines(self):
        topology = self.make_topology()
        invalid_face = LayoutFace(1, ((1, 1), (2, 9), (3, 3)))
        invalid = LayoutTopology(
            topology.points, topology.lines, (invalid_face,), topology.blocks
        )
        with self.assertRaisesRegex(ValueError, "unknown line"):
            invalid.validate()

    def test_face_sides_must_form_a_closed_loop(self):
        topology = self.make_topology()
        disconnected_lines = (
            topology.lines[0],
            LayoutLine(2, 1, 3, 3, ("auto", 0.1), 2),
            topology.lines[2],
        )
        invalid = LayoutTopology(
            topology.points,
            disconnected_lines,
            topology.faces,
            topology.blocks,
        )
        with self.assertRaisesRegex(ValueError, "share exactly one endpoint"):
            invalid.validate()

    def test_block_must_reference_existing_faces(self):
        topology = self.make_topology()
        invalid = LayoutTopology(
            topology.points,
            topology.lines,
            topology.faces,
            (LayoutBlock(1, ((1, 9),)),),
        )
        with self.assertRaisesRegex(ValueError, "unknown face"):
            invalid.validate()

    def test_duplicate_entity_ids_are_rejected(self):
        topology = self.make_topology()
        invalid = LayoutTopology(
            topology.points + (topology.points[0],),
            topology.lines,
            topology.faces,
            topology.blocks,
        )
        with self.assertRaisesRegex(ValueError, "Duplicate point identifier"):
            invalid.validate()


if __name__ == "__main__":
    unittest.main()
