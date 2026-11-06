"""Unit tests for the MeshSnapshot interchange manifest contract."""

from __future__ import annotations

import copy
import unittest

from .snapshot import MeshSnapshotError, validate_snapshot


def sample_snapshot() -> dict:
    return {
        "schema": "oneflow.mesh-snapshot",
        "schemaVersion": 1,
        "dimension": 2,
        "nodes": [
            {"id": 10, "coordinates": [0.0, 0.0]},
            {"id": 20, "coordinates": [1.0, 0.0]},
            {"id": 30, "coordinates": [1.0, 1.0]},
            {"id": 40, "coordinates": [0.0, 1.0]},
        ],
        "cells": [
            {"id": 100, "type": "triangle3", "nodeIds": [10, 20, 30]},
            {"id": 200, "type": "triangle3", "nodeIds": [10, 30, 40]},
        ],
        "faces": [
            {"id": 1, "nodeIds": [10, 20], "ownerCellId": 100, "neighbourCellId": None},
            {"id": 2, "nodeIds": [20, 30], "ownerCellId": 100, "neighbourCellId": None},
            {"id": 3, "nodeIds": [30, 10], "ownerCellId": 100, "neighbourCellId": 200},
            {"id": 4, "nodeIds": [30, 40], "ownerCellId": 200, "neighbourCellId": None},
            {"id": 5, "nodeIds": [40, 10], "ownerCellId": 200, "neighbourCellId": None},
        ],
        "groups": [
            {"id": "fluid", "name": "Fluid", "kind": "region", "memberKind": "cell", "dimension": 2, "memberIds": [100, 200]},
            {"id": "wall", "name": "Wall", "kind": "boundary", "memberKind": "face", "dimension": 1, "memberIds": [1, 2, 4, 5]},
        ],
        "provenance": {"sourceFormat": "gmsh-msh-4.1-ascii"},
    }


class MeshSnapshotContractTest(unittest.TestCase):
    def test_accepts_sparse_source_node_and_cell_ids(self) -> None:
        validate_snapshot(sample_snapshot())

    def test_rejects_duplicate_node_ids(self) -> None:
        snapshot = sample_snapshot()
        snapshot["nodes"][1]["id"] = 10
        with self.assertRaisesRegex(MeshSnapshotError, "duplicate node id"):
            validate_snapshot(snapshot)

    def test_rejects_unknown_connectivity_reference(self) -> None:
        snapshot = sample_snapshot()
        snapshot["cells"][0]["nodeIds"][0] = 999
        with self.assertRaisesRegex(MeshSnapshotError, "unknown node ID"):
            validate_snapshot(snapshot)

    def test_rejects_wrong_coordinate_dimension(self) -> None:
        snapshot = sample_snapshot()
        snapshot["nodes"][0]["coordinates"].append(0.0)
        with self.assertRaisesRegex(MeshSnapshotError, "coordinates must contain 2"):
            validate_snapshot(snapshot)

    def test_rejects_non_finite_coordinates(self) -> None:
        snapshot = sample_snapshot()
        snapshot["nodes"][0]["coordinates"][0] = float("nan")
        with self.assertRaisesRegex(MeshSnapshotError, "finite numbers"):
            validate_snapshot(snapshot)

    def test_rejects_invalid_cell_arity(self) -> None:
        snapshot = sample_snapshot()
        snapshot["cells"][0]["nodeIds"].pop()
        with self.assertRaisesRegex(MeshSnapshotError, "invalid arity"):
            validate_snapshot(snapshot)

    def test_rejects_group_members_not_present_in_cells(self) -> None:
        snapshot = sample_snapshot()
        snapshot["groups"][0]["memberIds"].append(999)
        with self.assertRaisesRegex(MeshSnapshotError, "unknown cell ID"):
            validate_snapshot(snapshot)

    def test_rejects_face_with_unknown_owner_cell(self) -> None:
        snapshot = sample_snapshot()
        snapshot["faces"][0]["ownerCellId"] = 999
        with self.assertRaisesRegex(MeshSnapshotError, "ownerCellId"):
            validate_snapshot(snapshot)

    def test_boundary_group_references_faces_not_cells(self) -> None:
        snapshot = sample_snapshot()
        snapshot["groups"][1]["memberIds"] = [100]
        with self.assertRaisesRegex(MeshSnapshotError, "unknown face ID"):
            validate_snapshot(snapshot)

    def test_rejects_duplicate_group_ids(self) -> None:
        snapshot = sample_snapshot()
        snapshot["groups"].append(copy.deepcopy(snapshot["groups"][0]))
        with self.assertRaisesRegex(MeshSnapshotError, "duplicate group id"):
            validate_snapshot(snapshot)


if __name__ == "__main__":
    unittest.main()
