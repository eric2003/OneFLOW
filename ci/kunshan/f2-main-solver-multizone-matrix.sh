#!/bin/bash
# Multi-zone HIP scalability matrix.
#
# Runs the same partitioned case through the paired-timing runner for several
# (rank, DCU) tuples.  Resource bookkeeping follows the Kunshan DCU node shape:
# 32 cores and 4 DCUs per node, i.e. one DCU slot is paired with 8 cores
# (`max_dcus_per_node=4`, `max_cpus_per_node=32`; verify with
# `scnet.py queues --details`).  A leg is "spec-compliant" only when
#   cores_per_dcu == 8, i.e. cores == 8 * dcus,
# and the CPU baseline must come from the CPU-only partition (`kshcnormal`),
# not from an allocation that also holds a DCU.  Set
# `ONEFLOW_MZ_MATRIX_CPU_BASELINE=.../timing.tsv` to compare against a baseline
# measured that way; otherwise the in-allocation CPU leg is used and labelled
# as such.
set -euo pipefail

: "${ONEFLOW_MZ_MATRIX_SOURCE_CASE:?Set ONEFLOW_MZ_MATRIX_SOURCE_CASE}"
: "${ONEFLOW_MZ_MATRIX_BINARY:?Set ONEFLOW_MZ_MATRIX_BINARY}"
: "${ONEFLOW_MZ_MATRIX_RUN_ROOT:?Set ONEFLOW_MZ_MATRIX_RUN_ROOT}"
: "${ONEFLOW_MZ_MATRIX_TIMING_RUNNER:?Set ONEFLOW_MZ_MATRIX_TIMING_RUNNER}"
: "${ONEFLOW_MZ_MATRIX_RANK_VERIFIER:?Set ONEFLOW_MZ_MATRIX_RANK_VERIFIER}"

steps="${ONEFLOW_MZ_MATRIX_STEPS:-10}"
warmup="${ONEFLOW_MZ_MATRIX_WARMUP:-1}"
repeats="${ONEFLOW_MZ_MATRIX_REPEATS:-3}"
cfl="${ONEFLOW_MZ_MATRIX_CFL:-0.01}"
timeout_seconds="${ONEFLOW_MZ_MATRIX_TIMEOUT_SECONDS:-1200}"
ranks_list="${ONEFLOW_MZ_MATRIX_RANKS:-8 4}"
dcu_list="${ONEFLOW_MZ_MATRIX_DCUS:-1 4}"
fixture_metadata="${ONEFLOW_MZ_MATRIX_FIXTURE_METADATA:-${ONEFLOW_MZ_FIXTURE_METADATA:-}}"
# Optional extra HIP stage opt-ins for the HIP leg, for example
# "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE=1 ONEFLOW_ENABLE_UNS_HIP_BOUNDARY=1".
hip_stage_env="${ONEFLOW_MZ_MATRIX_HIP_STAGE_ENV:-}"
cores_per_task="${ONEFLOW_MZ_MATRIX_CORES_PER_TASK:-1}"
external_cpu_baseline="${ONEFLOW_MZ_MATRIX_CPU_BASELINE:-}"
allow_nonstandard="${ONEFLOW_MZ_MATRIX_ALLOW_NONSTANDARD:-0}"

test -d "$ONEFLOW_MZ_MATRIX_SOURCE_CASE"
test -x "$ONEFLOW_MZ_MATRIX_BINARY"
test -r "$ONEFLOW_MZ_MATRIX_TIMING_RUNNER"
test -r "$ONEFLOW_MZ_MATRIX_RANK_VERIFIER"
if [ -e "$ONEFLOW_MZ_MATRIX_RUN_ROOT" ]; then
    echo "refusing to overwrite $ONEFLOW_MZ_MATRIX_RUN_ROOT" >&2
    exit 2
fi
mkdir -p "$ONEFLOW_MZ_MATRIX_RUN_ROOT"

{
    echo "case=$ONEFLOW_MZ_MATRIX_SOURCE_CASE"
    echo "binary=$ONEFLOW_MZ_MATRIX_BINARY"
    echo "ranks=$ranks_list"
    echo "dcus=$dcu_list"
    echo "steps=$steps"
    echo "warmup=$warmup"
    echo "repeats=$repeats"
    echo "cfl=$cfl"
    echo "hip_stage_env=${hip_stage_env:-default-policy}"
    echo "cores_per_task=$cores_per_task"
    echo "slot_rule=one_dcu_slot_equals_8_cores"
    echo "external_cpu_baseline=${external_cpu_baseline:-in-allocation}"
    echo "fixture_metadata=${fixture_metadata:-unset}"
    echo "nodes=${SLURM_JOB_NUM_NODES:-1}"
    echo "hostlist=${SLURM_JOB_NODELIST:-unknown}"
    echo "mpirun=$(command -v mpirun)"
    echo "mpi=$(mpirun --version | head -1)"
} > "$ONEFLOW_MZ_MATRIX_RUN_ROOT/environment.txt"

visible_devices() {
    local count="$1"
    local list=""
    local index
    for (( index = 0; index < count; ++index )); do
        list="${list:+$list,}$index"
    done
    printf '%s' "$list"
}

run_leg() {
    local mode="$1"
    local ranks="$2"
    local dcus="$3"
    local run_root="$ONEFLOW_MZ_MATRIX_RUN_ROOT/$mode-r${ranks}-d${dcus}"

    (
        export ONEFLOW_MZ_SOURCE_CASE="$ONEFLOW_MZ_MATRIX_SOURCE_CASE"
        export ONEFLOW_MZ_BINARY="$ONEFLOW_MZ_MATRIX_BINARY"
        export ONEFLOW_MZ_RUN_ROOT="$run_root"
        export ONEFLOW_MZ_MODE="$mode"
        export ONEFLOW_MZ_RANKS="$ranks"
        export ONEFLOW_MZ_STEPS="$steps"
        export ONEFLOW_MZ_WARMUP="$warmup"
        export ONEFLOW_MZ_REPEATS="$repeats"
        export ONEFLOW_MZ_CFL="$cfl"
        export ONEFLOW_MZ_TIMEOUT_SECONDS="$timeout_seconds"
        export ONEFLOW_MZ_RANK_VERIFIER="$ONEFLOW_MZ_MATRIX_RANK_VERIFIER"
        if [ -n "$fixture_metadata" ]; then
            export ONEFLOW_MZ_FIXTURE_METADATA="$fixture_metadata"
        fi
        if [ "$mode" = hip ]; then
            local cores=$(( ranks * cores_per_task ))
            local cores_per_dcu=$(( cores / dcus ))
            if [ "$cores_per_dcu" -ne 8 ] && [ "$allow_nonstandard" != 1 ]; then
                echo "non-standard slot: cores=$cores dcus=$dcus cores_per_dcu=$cores_per_dcu (expected 8); set ONEFLOW_MZ_MATRIX_ALLOW_NONSTANDARD=1 to record it as a card-scaling experiment" >&2
                exit 3
            fi
            local devices
            devices="$(visible_devices "$dcus")"
            export HIP_VISIBLE_DEVICES="$devices"
            export ROCR_VISIBLE_DEVICES="$devices"
            if [ -n "$hip_stage_env" ]; then
                # shellcheck disable=SC2086
                export $hip_stage_env
            fi
        fi
        bash "$ONEFLOW_MZ_MATRIX_TIMING_RUNNER"
    ) > "$ONEFLOW_MZ_MATRIX_RUN_ROOT/$mode-r${ranks}-d${dcus}.log" 2>&1
    printf '%s\n' "$run_root/timing.tsv"
}

timing_files=()
for ranks in $ranks_list; do
    # One CPU baseline per rank count: adding DCUs must not change the CPU leg.
    timing_files+=("$(run_leg cpu "$ranks" 0)")
    for dcus in $dcu_list; do
        if [ "$dcus" -gt 0 ]; then
            timing_files+=("$(run_leg hip "$ranks" "$dcus")")
        fi
    done
done

python3 - "$ONEFLOW_MZ_MATRIX_RUN_ROOT" "${timing_files[@]}" <<'PY'
import csv
import os
import re
import statistics
import sys

root = sys.argv[1]
paths = sys.argv[2:]
rows = []
for path in paths:
    if not os.path.isfile(path):
        raise SystemExit(f"missing timing file: {path}")
    case = os.path.basename(os.path.dirname(path))
    with open(path, newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            if record.get("phase") != "formal":
                continue
            rows.append(
                {
                    "case": case,
                    "mode": record["mode"],
                    "ranks": int(record["ranks"]),
                    "elapsed_ms": float(record["elapsed_ms"]),
                    "exit_code": int(record["exit_code"]),
                }
            )

expected = {}
for path in paths:
    case = os.path.basename(os.path.dirname(path))
    with open(path, newline="") as handle:
        formal = [
            record
            for record in csv.DictReader(handle, delimiter="\t")
            if record.get("phase") == "formal"
        ]
    expected[case] = len(formal)

aggregate = {}
for row in rows:
    aggregate.setdefault((row["case"], row["mode"], row["ranks"]), []).append(row)

lines = [
    "case\tmode\tranks\tdcus\tranks_per_dcu\tcores\tcores_per_dcu\trepeats\tmean_ms\tmin_ms\tmax_ms\texit_code"
]
summary = {}
for (case, mode, ranks), values in sorted(aggregate.items()):
    dcus = 0
    match = re.search(r"-d(\d+)$", case)
    if match:
        dcus = int(match.group(1))
    if any(value["exit_code"] != 0 for value in values):
        raise SystemExit(f"{case}: non-zero workload exit code")
    if len(values) != expected[case]:
        raise SystemExit(
            f"{case}: aggregated {len(values)} repeats, expected {expected[case]}"
        )
    elapsed = [value["elapsed_ms"] for value in values]
    mean = statistics.fmean(elapsed)
    summary[(mode, ranks, dcus)] = mean
    cores = ranks
    lines.append(
        "\t".join(
            [
                case,
                mode,
                str(ranks),
                str(dcus),
                "-" if mode == "cpu" else f"{ranks / dcus:g}",
                str(cores),
                "-" if mode == "cpu" else f"{cores / dcus:g}",
                str(len(elapsed)),
                f"{mean:.6f}",
                f"{min(elapsed):.6f}",
                f"{max(elapsed):.6f}",
                "0",
            ]
        )
    )

ratio_lines = [
    "ranks\tdcus\ttranks_per_dcu\tcores_per_dcu\tcpu_mean_ms\thip_mean_ms\tcpu_over_hip"
]
rank_values = sorted({key[1] for key in summary if key[0] == "cpu"})
dcu_values = sorted({key[2] for key in summary if key[0] == "hip" and key[2] > 0})
for ranks in rank_values:
    cpu_mean = summary.get(("cpu", ranks, 0))
    if cpu_mean is None:
        continue
    for dcus in dcu_values:
        hip_mean = summary.get(("hip", ranks, dcus))
        if hip_mean is None:
            continue
        ratio_lines.append(
            "\t".join(
                [
                    str(ranks),
                    str(dcus),
                    f"{ranks / dcus:g}",
                    f"{ranks / dcus:g}",
                    f"{cpu_mean:.6f}",
                    f"{hip_mean:.6f}",
                    f"{cpu_mean / hip_mean:.6f}",
                ]
            )
        )

with open(os.path.join(root, "matrix.tsv"), "w") as handle:
    handle.write("\n".join(lines) + "\n")
with open(os.path.join(root, "matrix-ratio.tsv"), "w") as handle:
    handle.write("\n".join(ratio_lines) + "\n")
print("\n".join(ratio_lines))
PY

printf 'MULTIZONE_MATRIX_PASS root=%s\n' "$ONEFLOW_MZ_MATRIX_RUN_ROOT" | tee "$ONEFLOW_MZ_MATRIX_RUN_ROOT/result.txt"
