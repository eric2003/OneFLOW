#!/usr/bin/env python3
"""Merge runtime interface shards into the v1 multi-zone fixture manifest.

The solver writes one shard per MPI rank after resolving its actual
InterFace/FaceTopo endpoint maps.  This tool deliberately refuses to infer
topology from aggregate face counts: every send endpoint must have a matching
receive endpoint on the neighboring zone and vice versa.
"""

import argparse
import json
import pathlib
import sys
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Tuple


SHARD_SCHEMA = "oneflow.multizone.fixture/shard/v1"
MANIFEST_SCHEMA = "oneflow.multizone.fixture/v1"
MEMORY_SPACES = {"HostStaged", "DeviceDirect"}


def _integer(value: Any, field: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{field} must be an integer >= {minimum}")
    return value


def _load_shard(path: pathlib.Path) -> Dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: cannot read JSON shard: {error}") from error
    if not isinstance(value, dict) or value.get("schema") != SHARD_SCHEMA:
        raise ValueError(f"{path}: schema must be {SHARD_SCHEMA!r}")
    _integer(value.get("rank"), f"{path}.rank")
    ranks = _integer(value.get("ranks"), f"{path}.ranks", 1)
    if value["rank"] >= ranks:
        raise ValueError(f"{path}: rank is outside shard rank count")
    if not isinstance(value.get("zones"), list):
        raise ValueError(f"{path}.zones must be an array")
    if not isinstance(value.get("endpoints"), list):
        raise ValueError(f"{path}.endpoints must be an array")
    return value


def _endpoint_key(item: Dict[str, Any]) -> Tuple[int, int, int, int]:
    return (
        item["rank"],
        item["zone"],
        item["neighbor_rank"],
        item["neighbor_zone"],
    )


def merge_shards(shards: Iterable[Dict[str, Any]]) -> Dict[str, Any]:
    loaded = list(shards)
    if not loaded:
        raise ValueError("at least one shard is required")
    expected_ranks = loaded[0]["ranks"]
    if any(shard["ranks"] != expected_ranks for shard in loaded):
        raise ValueError("all shards must use the same ranks value")

    zones: Dict[int, Dict[str, int]] = {}
    seen_ranks = set()
    sends: Dict[Tuple[int, int, int, int, int], Dict[str, Any]] = {}
    receives: Dict[Tuple[int, int, int, int, int], Dict[str, Any]] = {}

    for shard in loaded:
        rank = shard["rank"]
        if rank in seen_ranks:
            raise ValueError(f"duplicate shard rank {rank}")
        seen_ranks.add(rank)
        for index, raw_zone in enumerate(shard["zones"]):
            if not isinstance(raw_zone, dict):
                raise ValueError(f"zones[{index}] must be an object")
            zone = _integer(raw_zone.get("zone"), f"zones[{index}].zone")
            owner = _integer(raw_zone.get("rank"), f"zones[{index}].rank")
            owned = _integer(
                raw_zone.get("owned_cells"), f"zones[{index}].owned_cells", 1
            )
            ghost = _integer(
                raw_zone.get("ghost_cells"), f"zones[{index}].ghost_cells"
            )
            if owner != rank:
                raise ValueError(f"zone {zone} is not owned by its shard rank")
            candidate = {
                "zone": zone,
                "rank": owner,
                "owned_cells": owned,
                "ghost_cells": ghost,
            }
            if zone in zones and zones[zone] != candidate:
                raise ValueError(f"zone {zone} has conflicting ownership metadata")
            zones[zone] = candidate

        for index, raw_endpoint in enumerate(shard["endpoints"]):
            if not isinstance(raw_endpoint, dict):
                raise ValueError(f"endpoints[{index}] must be an object")
            endpoint = dict(raw_endpoint)
            for field in (
                "zone",
                "rank",
                "neighbor_zone",
                "neighbor_rank",
                "ordinal",
                "face",
                "cell",
                "equations",
            ):
                _integer(endpoint.get(field), f"endpoints[{index}].{field}")
            if endpoint["rank"] != rank:
                raise ValueError(f"endpoint {index} rank disagrees with shard")
            if endpoint.get("direction") not in {"send", "receive"}:
                raise ValueError(f"endpoints[{index}].direction is invalid")
            if endpoint.get("role") not in {"OwnedCell", "GhostCell"}:
                raise ValueError(f"endpoints[{index}].role is invalid")
            if endpoint["memory_space"] not in MEMORY_SPACES:
                raise ValueError(f"endpoints[{index}].memory_space is invalid")
            key = (
                endpoint["zone"],
                endpoint["neighbor_zone"],
                endpoint["ordinal"],
                endpoint["rank"],
                endpoint["neighbor_rank"],
            )
            target = sends if endpoint["direction"] == "send" else receives
            if key in target:
                raise ValueError(f"duplicate endpoint key {key}")
            target[key] = endpoint

    if seen_ranks != set(range(expected_ranks)):
        raise ValueError(
            f"shards must cover ranks 0..{expected_ranks - 1}; "
            f"missing={sorted(set(range(expected_ranks)) - seen_ranks)}"
        )
    if set(zones) != set(range(expected_ranks)):
        raise ValueError(
            "fixture must contain one zone per rank; "
            f"zones={sorted(zones)} ranks={expected_ranks}"
        )

    interfaces: List[Dict[str, Any]] = []
    for send_key, send in sorted(sends.items()):
        zone, neighbor_zone, ordinal, rank, neighbor_rank = send_key
        receive_key = (
            neighbor_zone,
            zone,
            ordinal,
            neighbor_rank,
            rank,
        )
        receive = receives.get(receive_key)
        if receive is None:
            raise ValueError(f"send endpoint {send_key} has no reciprocal receive")
        reverse_send = sends.get(
            (neighbor_zone, zone, ordinal, neighbor_rank, rank)
        )
        reverse_receive = receives.get(
            (zone, neighbor_zone, ordinal, rank, neighbor_rank)
        )
        if reverse_send is None or reverse_receive is None:
            raise ValueError(f"endpoint pair {send_key} is not bidirectional")
        if reverse_send["face"] != receive["face"]:
            raise ValueError(
                f"endpoint pair {send_key} disagrees on target face ID"
            )
        if reverse_receive["face"] != send["face"]:
            raise ValueError(
                f"endpoint pair {send_key} disagrees on source face ID"
            )
        if send["role"] != "OwnedCell" or receive["role"] != "GhostCell":
            raise ValueError(f"endpoint pair {send_key} has invalid roles")
        if send["equations"] != receive["equations"]:
            raise ValueError(f"endpoint pair {send_key} equation counts differ")
        if send["memory_space"] != receive["memory_space"]:
            raise ValueError(f"endpoint pair {send_key} memory spaces differ")
        equations = send["equations"]
        interfaces.append(
            {
                "source": {
                    "rank": rank,
                    "zone": zone,
                    "face": send["face"],
                    "cell": send["cell"],
                    "role": "OwnedCell",
                },
                "target": {
                    "rank": neighbor_rank,
                    "zone": neighbor_zone,
                    "face": receive["face"],
                    "cell": receive["cell"],
                    "role": "GhostCell",
                },
                "equations": equations,
                "send_values": equations,
                "receive_values": equations,
                "memory_space": send["memory_space"],
            }
        )

    if len(sends) != len(receives):
        raise ValueError(
            f"send/receive endpoint counts differ: {len(sends)} vs {len(receives)}"
        )
    return {
        "schema": MANIFEST_SCHEMA,
        "ranks": expected_ranks,
        "zones": [zones[zone] for zone in sorted(zones)],
        "interfaces": interfaces,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        paths = sorted(args.shards.glob("rank-*.json"))
        manifest = merge_shards(_load_shard(path) for path in paths)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    except (OSError, ValueError) as error:
        print(f"MULTIZONE_FIXTURE_MERGE_FAIL error={error}", file=sys.stderr)
        return 1
    print(
        "MULTIZONE_FIXTURE_MERGE_PASS "
        f"output={args.output} ranks={manifest['ranks']} "
        f"zones={len(manifest['zones'])} interfaces={len(manifest['interfaces'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
