#!/usr/bin/env python3
"""Tests for the external multi-zone fixture metadata gate."""

import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).with_name("verify-multizone-fixture.py")
SPEC = importlib.util.spec_from_file_location("verify_multizone_fixture", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
MERGE_SCRIPT = pathlib.Path(__file__).with_name("merge-multizone-fixture.py")
MERGE_SPEC = importlib.util.spec_from_file_location(
    "merge_multizone_fixture", MERGE_SCRIPT
)
assert MERGE_SPEC and MERGE_SPEC.loader
MERGE_MODULE = importlib.util.module_from_spec(MERGE_SPEC)
MERGE_SPEC.loader.exec_module(MERGE_MODULE)


def manifest():
    return {
        "schema": "oneflow.multizone.fixture/v1",
        "ranks": 2,
        "zones": [
            {"zone": 0, "rank": 0, "owned_cells": 10, "ghost_cells": 2},
            {"zone": 1, "rank": 1, "owned_cells": 12, "ghost_cells": 2},
        ],
        "interfaces": [
            {
                "source": {
                    "rank": 0,
                    "zone": 0,
                    "face": 3,
                    "cell": 4,
                    "role": "OwnedCell",
                },
                "target": {
                    "rank": 1,
                    "zone": 1,
                    "face": 5,
                    "cell": 12,
                    "role": "GhostCell",
                },
                "equations": 5,
                "send_values": 5,
                "receive_values": 5,
                "memory_space": "HostStaged",
            },
            {
                "source": {
                    "rank": 1,
                    "zone": 1,
                    "face": 5,
                    "cell": 6,
                    "role": "OwnedCell",
                },
                "target": {
                    "rank": 0,
                    "zone": 0,
                    "face": 3,
                    "cell": 10,
                    "role": "GhostCell",
                },
                "equations": 5,
                "send_values": 5,
                "receive_values": 5,
                "memory_space": "HostStaged",
            },
        ],
    }


class MultiZoneFixtureVerifierTest(unittest.TestCase):
    def test_valid_manifest(self):
        summary = MODULE.validate_manifest(manifest(), 2)
        self.assertEqual(summary, {"ranks": 2, "zones": 2, "interfaces": 2})

    def test_rejects_self_neighbor(self):
        value = manifest()
        value["interfaces"][0]["target"]["rank"] = 0
        with self.assertRaisesRegex(ValueError, "cannot exchange with itself"):
            MODULE.validate_manifest(value, 2)

    def test_rejects_wrong_endpoint_role(self):
        value = manifest()
        value["interfaces"][0]["source"]["role"] = "GhostCell"
        with self.assertRaisesRegex(ValueError, "OwnedCell source"):
            MODULE.validate_manifest(value, 2)

    def test_allows_zero_length_message(self):
        value = manifest()
        for interface in value["interfaces"]:
            interface["send_values"] = 0
            interface["receive_values"] = 0
        summary = MODULE.validate_manifest(value, 2)
        self.assertEqual(summary["interfaces"], 2)

    def test_rejects_target_cell_outside_zone(self):
        value = manifest()
        value["interfaces"][0]["target"]["cell"] = 14
        with self.assertRaisesRegex(ValueError, "outside target zone"):
            MODULE.validate_manifest(value, 2)

    def test_cli_reports_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = pathlib.Path(temporary) / "fixture.json"
            path.write_text(json.dumps(manifest()), encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--metadata",
                    str(path),
                    "--expected-ranks",
                    "2",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("MULTIZONE_FIXTURE_PASS", result.stdout)

    def test_runtime_shards_merge_without_inference(self):
        send = {
            "schema": "oneflow.multizone.fixture/shard/v1",
            "rank": 0,
            "ranks": 2,
            "zones": [{"zone": 0, "rank": 0, "owned_cells": 10, "ghost_cells": 2}],
            "endpoints": [
                {
                    "zone": 0,
                    "rank": 0,
                    "neighbor_zone": 1,
                    "neighbor_rank": 1,
                    "ordinal": 0,
                    "direction": "send",
                    "face": 3,
                    "cell": 4,
                    "role": "OwnedCell",
                    "equations": 5,
                    "memory_space": "HostStaged",
                },
                {
                    "zone": 0,
                    "rank": 0,
                    "neighbor_zone": 1,
                    "neighbor_rank": 1,
                    "ordinal": 0,
                    "direction": "receive",
                    "face": 3,
                    "cell": 10,
                    "role": "GhostCell",
                    "equations": 5,
                    "memory_space": "HostStaged",
                },
            ],
        }
        receive = {
            "schema": "oneflow.multizone.fixture/shard/v1",
            "rank": 1,
            "ranks": 2,
            "zones": [{"zone": 1, "rank": 1, "owned_cells": 12, "ghost_cells": 2}],
            "endpoints": [
                {
                    "zone": 1,
                    "rank": 1,
                    "neighbor_zone": 0,
                    "neighbor_rank": 0,
                    "ordinal": 0,
                    "direction": "receive",
                    "face": 5,
                    "cell": 12,
                    "role": "GhostCell",
                    "equations": 5,
                    "memory_space": "HostStaged",
                },
                {
                    "zone": 1,
                    "rank": 1,
                    "neighbor_zone": 0,
                    "neighbor_rank": 0,
                    "ordinal": 0,
                    "direction": "send",
                    "face": 5,
                    "cell": 6,
                    "role": "OwnedCell",
                    "equations": 5,
                    "memory_space": "HostStaged",
                },
            ],
        }
        merged = MERGE_MODULE.merge_shards([send, receive])
        self.assertEqual(merged["schema"], "oneflow.multizone.fixture/v1")
        self.assertEqual(
            MODULE.validate_manifest(merged, 2),
            {"ranks": 2, "zones": 2, "interfaces": 2},
        )

    def test_runtime_shards_reject_missing_reciprocal(self):
        value = {
            "schema": "oneflow.multizone.fixture/shard/v1",
            "rank": 0,
            "ranks": 1,
            "zones": [{"zone": 0, "rank": 0, "owned_cells": 2, "ghost_cells": 1}],
            "endpoints": [
                {
                    "zone": 0,
                    "rank": 0,
                    "neighbor_zone": 0,
                    "neighbor_rank": 0,
                    "ordinal": 0,
                    "direction": "send",
                    "face": 0,
                    "cell": 0,
                    "role": "OwnedCell",
                    "equations": 1,
                    "memory_space": "HostStaged",
                }
            ],
        }
        with self.assertRaisesRegex(ValueError, "no reciprocal receive"):
            MERGE_MODULE.merge_shards([value])


if __name__ == "__main__":
    unittest.main()
