"""Unit tests for the MeshSnapshot interchange manifest contract."""

from __future__ import annotations

import copy
import unittest

from mesh_contracts.snapshot import MeshSnapshotError, validate_snapshot


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
        "groups": [
            {"id": "fluid", "name": "Fluid", "kind": "region", "dimension": 2, "memberIds": [100, 200]},
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

    def test_rejects_duplicate_group_ids(self) -> None:
        snapshot = sample_snapshot()
        snapshot["groups"].append(copy.deepcopy(snapshot["groups"][0]))
        with self.assertRaisesRegex(MeshSnapshotError, "duplicate group id"):
            validate_snapshot(snapshot)


if __name__ == "__main__":
    unittest.main()
