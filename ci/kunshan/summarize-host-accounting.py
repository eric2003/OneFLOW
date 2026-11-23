#!/usr/bin/env python3
"""Validate MG wall accounting and report startup-corrected per-rank scopes.

Exclusive times partition MG::Run, including device/MPI waits. They are not
CPU execution times. MPI rank times must never be summed into wall time.
"""

import argparse
import csv
import glob
import json
import math
import pathlib
import re
import sys


def read_profiles(pattern, ranks, steps):
    profiles = {}
    for filename in sorted(glob.glob(pattern)):
        match = re.search(r"-(\d+)\.tsv\.host\.tsv$", filename)
        if not match:
            raise ValueError("invalid rank filename: " + filename)
        rank = int(match.group(1))
        if rank in profiles:
            raise ValueError("duplicate rank: " + str(rank))
        rows = {}
        with open(filename, newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames != ["scope", "count", "inclusive_ms", "exclusive_ms"]:
                raise ValueError("unexpected host accounting header: " + filename)
            for row in reader:
                name = row["scope"]
                count = int(row["count"])
                inclusive = float(row["inclusive_ms"])
                exclusive = float(row["exclusive_ms"])
                if (not name or name in rows or count <= 0 or
                        not all(math.isfinite(v) and v >= 0 for v in (inclusive, exclusive)) or
                        exclusive > inclusive + 1e-8):
                    raise ValueError("invalid/duplicate accounting row: " + filename)
                rows[name] = dict(count=count, inclusive_ms=inclusive, exclusive_ms=exclusive)
        if rows.get("mg_run", {}).get("count") != 1:
            raise ValueError("expected exactly one mg_run root: " + filename)
        if rows.get("mg_inner_iteration", {}).get("count") != steps:
            raise ValueError("unexpected inner iteration count: " + filename)
        root = rows["mg_run"]["inclusive_ms"]
        if root <= 0 or abs(sum(r["exclusive_ms"] for r in rows.values()) - root) > max(1e-6, root * 1e-9):
            raise ValueError("exclusive scope sum does not close to mg_run: " + filename)
        profiles[rank] = rows
    if set(profiles) != set(range(ranks)):
        raise ValueError("missing or extra MPI rank files")
    return profiles


def summarize(long_profiles, short_profiles, steps, short_steps, wall_ms, short_wall_ms):
    delta_steps = steps - short_steps
    if delta_steps <= 0 or not all(math.isfinite(v) and v > 0 for v in (wall_ms, short_wall_ms)):
        raise ValueError("invalid measurement basis")
    wall_per_step = (wall_ms - short_wall_ms) / delta_steps
    if wall_per_step <= 0:
        raise ValueError("non-positive startup-corrected wall time")
    output = []
    for rank, long_rows in sorted(long_profiles.items()):
        short_rows = short_profiles[rank]
        names = set(long_rows) | set(short_rows)
        scopes = []
        for name in names:
            total = long_rows.get(name, {}).get("exclusive_ms", 0.0)
            short = short_rows.get(name, {}).get("exclusive_ms", 0.0)
            # Small signed differences can occur for one-time work/noise.
            scopes.append(dict(scope=name, exclusive_ms_per_step=(total - short) / delta_steps))
        root_per_step = (long_rows["mg_run"]["inclusive_ms"] -
                         short_rows["mg_run"]["inclusive_ms"]) / delta_steps
        if root_per_step <= 0:
            raise ValueError("non-positive startup-corrected MG time")
        output.append(dict(rank=rank, mg_ms_per_step=root_per_step,
                           scopes=sorted(scopes, key=lambda r: -r["exclusive_ms_per_step"])))
    # The slowest rank is a wall-time diagnostic; do not sum rank durations.
    critical = max(output, key=lambda r: r["mg_ms_per_step"])
    gap = abs(wall_per_step - critical["mg_ms_per_step"]) / wall_per_step
    return dict(status="PASS" if gap < 0.10 else "FAIL", wall_ms_per_step=wall_per_step,
                critical_rank=critical["rank"], unexplained_fraction=gap,
                steps=steps, short_steps=short_steps, ranks=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern", required=True)
    parser.add_argument("--short-pattern", required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--short-steps", type=int, default=1)
    parser.add_argument("--wall-ms", type=float, required=True)
    parser.add_argument("--short-wall-ms", type=float, required=True)
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()
    try:
        if args.ranks <= 0 or args.short_steps <= 0 or args.steps <= args.short_steps:
            raise ValueError("invalid ranks/step counts")
        result = summarize(read_profiles(args.pattern, args.ranks, args.steps),
                           read_profiles(args.short_pattern, args.ranks, args.short_steps),
                           args.steps, args.short_steps, args.wall_ms, args.short_wall_ms)
        content = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output:
            args.output.write_text(content, encoding="utf-8")
        print(content, end="")
        return 0 if result["status"] == "PASS" else 1
    except (ValueError, OSError, KeyError, TypeError) as error:
        print("HOST_ACCOUNTING_FAIL: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
