#!/usr/bin/env python3
"""Validate the small topology manifest paired with an external multi-zone case.

The mesh itself remains outside Git.  This manifest is intentionally limited to
ownership and exchange metadata so runners can reject a single-zone case or an
incomplete endpoint map before spending Slurm time on a solver invocation.
"""

import argparse
import json
import pathlib
import sys
from typing import Any, Dict, Tuple


SCHEMA = "oneflow.multizone.fixture/v1"
MEMORY_SPACES = {"HostStaged", "DeviceDirect"}


def _positive_int(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{field} must be a positive integer")
    return value


def _nonnegative_int(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field} must be a non-negative integer")
    return value


def _dict(value: Any, field: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{field} must be an object")
    return value


def _required(obj: Dict[str, Any], field: str) -> Any:
    if field not in obj:
        raise ValueError(f"missing {field}")
    return obj[field]


def _validate_zones(
    zones: Any, expected_ranks: int
) -> Tuple[Dict[int, Dict[str, Any]], Dict[int, int]]:
    if not isinstance(zones, list) or not zones:
        raise ValueError("zones must be a non-empty array")
    by_rank: Dict[int, Dict[str, Any]] = {}
    by_zone: Dict[int, int] = {}
    for index, raw_zone in enumerate(zones):
        zone = _dict(raw_zone, f"zones[{index}]")
        zone_id = _nonnegative_int(_required(zone, "zone"), f"zones[{index}].zone")
        rank = _nonnegative_int(_required(zone, "rank"), f"zones[{index}].rank")
        owned = _positive_int(
            _required(zone, "owned_cells"), f"zones[{index}].owned_cells"
        )
        ghost = _nonnegative_int(
            _required(zone, "ghost_cells"), f"zones[{index}].ghost_cells"
        )
        if rank in by_rank:
            raise ValueError(f"duplicate zone owner rank {rank}")
        if zone_id in by_zone:
            raise ValueError(f"duplicate zone id {zone_id}")
        by_rank[rank] = {"zone": zone_id, "owned_cells": owned, "ghost_cells": ghost}
        by_zone[zone_id] = rank

    expected = set(range(expected_ranks))
    actual = set(by_rank)
    if actual != expected:
        raise ValueError(
            f"zone owners must be ranks 0..{expected_ranks - 1}; "
            f"missing={sorted(expected - actual)} extra={sorted(actual - expected)}"
        )
    return by_rank, by_zone


def _endpoint(raw: Any, field: str) -> Dict[str, Any]:
    value = _dict(raw, field)
    role = _required(value, "role")
    if role not in {"OwnedCell", "GhostCell"}:
        raise ValueError(f"{field}.role must be OwnedCell or GhostCell")
    for name in ("rank", "zone", "face", "cell"):
        _nonnegative_int(_required(value, name), f"{field}.{name}")
    return value


def _validate_interfaces(
    interfaces: Any,
    zones_by_rank: Dict[int, Dict[str, Any]],
) -> int:
    if not isinstance(interfaces, list):
        raise ValueError("interfaces must be an array")
    pairs: Dict[Tuple[int, int], int] = {}
    for index, raw_interface in enumerate(interfaces):
        interface = _dict(raw_interface, f"interfaces[{index}]")
        source = _endpoint(_required(interface, "source"), f"interfaces[{index}].source")
        target = _endpoint(_required(interface, "target"), f"interfaces[{index}].target")
        if source["role"] != "OwnedCell" or target["role"] != "GhostCell":
            raise ValueError(
                f"interfaces[{index}] must map OwnedCell source to GhostCell target"
            )
        source_rank, target_rank = source["rank"], target["rank"]
        if source_rank == target_rank:
            raise ValueError(f"interfaces[{index}] cannot exchange with itself")
        if source_rank not in zones_by_rank or target_rank not in zones_by_rank:
            raise ValueError(f"interfaces[{index}] references an unknown rank")
        if source["zone"] != zones_by_rank[source_rank]["zone"]:
            raise ValueError(f"interfaces[{index}].source.zone disagrees with zone owner")
        if target["zone"] != zones_by_rank[target_rank]["zone"]:
            raise ValueError(f"interfaces[{index}].target.zone disagrees with zone owner")
        if source["cell"] >= zones_by_rank[source_rank]["owned_cells"]:
            raise ValueError(f"interfaces[{index}].source.cell is not owned by source zone")
        if target["cell"] < zones_by_rank[target_rank]["owned_cells"]:
            raise ValueError(f"interfaces[{index}].target.cell is not a ghost cell")
        if target["cell"] >= (
            zones_by_rank[target_rank]["owned_cells"]
            + zones_by_rank[target_rank]["ghost_cells"]
        ):
            raise ValueError(f"interfaces[{index}].target.cell is outside target zone")

        equations = _positive_int(
            _required(interface, "equations"), f"interfaces[{index}].equations"
        )
        send_values = _nonnegative_int(
            _required(interface, "send_values"), f"interfaces[{index}].send_values"
        )
        receive_values = _nonnegative_int(
            _required(interface, "receive_values"),
            f"interfaces[{index}].receive_values",
        )
        if send_values != receive_values:
            raise ValueError(f"interfaces[{index}] send/receive value counts differ")
        if send_values % equations != 0:
            raise ValueError(
                f"interfaces[{index}].send_values is not equation-major aligned"
            )
        memory_space = _required(interface, "memory_space")
        if memory_space not in MEMORY_SPACES:
            raise ValueError(
                f"interfaces[{index}].memory_space must be HostStaged or DeviceDirect"
            )
        key = (source_rank, target_rank)
        pairs[key] = pairs.get(key, 0) + 1

    if not interfaces:
        raise ValueError("interfaces must contain at least one exchange endpoint")
    for (source_rank, target_rank), count in pairs.items():
        if pairs.get((target_rank, source_rank)) != count:
            raise ValueError(
                f"interface pair {source_rank}<->{target_rank} is not reciprocal"
            )
    return len(interfaces)


def validate_manifest(manifest: Any, expected_ranks: int) -> Dict[str, int]:
    root = _dict(manifest, "manifest")
    if root.get("schema") != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}")
    manifest_ranks = _positive_int(root.get("ranks"), "ranks")
    if manifest_ranks != expected_ranks:
        raise ValueError(
            f"manifest ranks={manifest_ranks} does not match expected ranks={expected_ranks}"
        )
    zones_by_rank, _ = _validate_zones(_required(root, "zones"), expected_ranks)
    interface_count = _validate_interfaces(
        _required(root, "interfaces"), zones_by_rank
    )
    return {
        "ranks": expected_ranks,
        "zones": len(zones_by_rank),
        "interfaces": interface_count,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=pathlib.Path)
    parser.add_argument("--expected-ranks", required=True, type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_ranks <= 0:
        print("expected ranks must be positive", file=sys.stderr)
        return 2
    try:
        with args.metadata.open(encoding="utf-8") as stream:
            summary = validate_manifest(json.load(stream), args.expected_ranks)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"MULTIZONE_FIXTURE_FAIL metadata={args.metadata} error={error}", file=sys.stderr)
        return 1
    print(
        "MULTIZONE_FIXTURE_PASS "
        f"metadata={args.metadata} ranks={summary['ranks']} "
        f"zones={summary['zones']} interfaces={summary['interfaces']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
