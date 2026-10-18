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
            points=(LayoutPoint(1, 0.0, 0.0, 0.0), LayoutPoint(2, 1.0, 0.0, 0.0)),
            lines=(LayoutLine(1, 1, 2, 3, ["auto", 0.1], 2),),
            faces=(LayoutFace(1, ((1, 1),)),),
            blocks=(LayoutBlock(1, ((1, 1),)),),
        )

    def test_valid_point_line_face_block_references(self):
        self.make_topology().validate()

    def test_line_must_reference_existing_points(self):
        topology = self.make_topology()
        invalid = LayoutTopology(topology.points, (LayoutLine(1, 1, 9, 3, [], 2),), topology.faces, topology.blocks)
        with self.assertRaisesRegex(ValueError, "unknown point"):
            invalid.validate()

    def test_face_must_reference_existing_lines(self):
        topology = self.make_topology()
        invalid = LayoutTopology(topology.points, topology.lines, (LayoutFace(1, ((1, 9),)),), topology.blocks)
        with self.assertRaisesRegex(ValueError, "unknown line"):
            invalid.validate()

    def test_block_must_reference_existing_faces(self):
        topology = self.make_topology()
        invalid = LayoutTopology(topology.points, topology.lines, topology.faces, (LayoutBlock(1, ((1, 9),)),))
        with self.assertRaisesRegex(ValueError, "unknown face"):
            invalid.validate()

    def test_duplicate_entity_ids_are_rejected(self):
        topology = self.make_topology()
        invalid = LayoutTopology(topology.points + (topology.points[0],), topology.lines, topology.faces, topology.blocks)
        with self.assertRaisesRegex(ValueError, "Duplicate point identifier"):
            invalid.validate()


if __name__ == "__main__":
    unittest.main()
