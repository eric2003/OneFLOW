#!/bin/bash
# Restart-equivalence gate for the production 3D solver.
#
# For each selected backend, compare
#   * one continuous run of N steps, with
#   * a run of S steps, a restart from the written flow file, and S more steps,
# using the same case, fixed CFL, RK scheme and backend.  The final flow file
# must match within machine precision, and the final residual/aero rows within
# the documented tolerances.
#
# Environment:
#   ONEFLOW_RESTART_SOURCE_CASE   prepared case directory (required)
#   ONEFLOW_RESTART_CPU_BINARY    legacy CPU binary (required)
#   ONEFLOW_RESTART_HIP_BINARY    HIP binary (required; may equal the CPU one)
#   ONEFLOW_RESTART_RUN_ROOT      append-only output directory (required)
#   ONEFLOW_RESTART_BACKENDS      "CPU HIP" (default)
#   ONEFLOW_RESTART_STEPS         total steps (default 100)
#   ONEFLOW_RESTART_SPLIT         restart step (default 50)
#   ONEFLOW_RESTART_CFL           fixed CFL for both ramp ends (default 0.01)
#   ONEFLOW_RESTART_TIMEOUT       per-run timeout in seconds (default 1800)
#   ONEFLOW_RESTART_FLOW_REL_TOL  flow relative tolerance (default 1e-12)
#   ONEFLOW_RESTART_RES_REL_TOL   residual relative tolerance (default 1e-12)
#   ONEFLOW_RESTART_AERO_REL_TOL  aero relative tolerance (default 1e-10)
#
# The script intentionally does not contain cluster paths, hostnames or job
# metadata.  The case directory owns the grid, boundary map and system tree.
# Scope: one rank / one local zone.  Multi-zone restart needs per-rank flow
# files and an MPI launch; adding that without a contract would hide failures.
set -euo pipefail

: "${ONEFLOW_RESTART_SOURCE_CASE:?Set ONEFLOW_RESTART_SOURCE_CASE}"
: "${ONEFLOW_RESTART_CPU_BINARY:?Set ONEFLOW_RESTART_CPU_BINARY}"
: "${ONEFLOW_RESTART_HIP_BINARY:?Set ONEFLOW_RESTART_HIP_BINARY}"
: "${ONEFLOW_RESTART_RUN_ROOT:?Set ONEFLOW_RESTART_RUN_ROOT}"

backends="${ONEFLOW_RESTART_BACKENDS:-CPU HIP}"
steps="${ONEFLOW_RESTART_STEPS:-100}"
split="${ONEFLOW_RESTART_SPLIT:-50}"
cfl="${ONEFLOW_RESTART_CFL:-0.01}"
timeout_seconds="${ONEFLOW_RESTART_TIMEOUT:-1800}"
flow_tolerance="${ONEFLOW_RESTART_FLOW_REL_TOL:-1e-12}"
residual_tolerance="${ONEFLOW_RESTART_RES_REL_TOL:-1e-12}"
aero_tolerance="${ONEFLOW_RESTART_AERO_REL_TOL:-1e-10}"
rank_verifier="${ONEFLOW_RESTART_RANK_VERIFIER:-$(dirname "$0")/verify-rank-utilization.py}"

for value_name in steps split timeout_seconds; do
    value="${!value_name}"
    case "$value" in
        ''|*[!0-9]*|0) echo "invalid $value_name=$value" >&2; exit 2 ;;
    esac
done
if [ "$split" -ge "$steps" ]; then
    echo "ONEFLOW_RESTART_SPLIT must be smaller than ONEFLOW_RESTART_STEPS" >&2
    exit 2
fi
if ! awk -v value="$cfl" 'BEGIN {
    valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/;
    exit !(valid && value + 0.0 > 0.0);
}'; then
    echo "invalid ONEFLOW_RESTART_CFL=$cfl" >&2
    exit 2
fi

test -d "$ONEFLOW_RESTART_SOURCE_CASE"
test -x "$ONEFLOW_RESTART_CPU_BINARY"
test -x "$ONEFLOW_RESTART_HIP_BINARY"
test -r "$rank_verifier"
if [ -e "$ONEFLOW_RESTART_RUN_ROOT" ]; then
    echo "refusing to overwrite $ONEFLOW_RESTART_RUN_ROOT" >&2
    exit 2
fi
mkdir -p "$ONEFLOW_RESTART_RUN_ROOT"

prepare_case() {
    local case_dir="$1"
    cp -a "$ONEFLOW_RESTART_SOURCE_CASE" "$case_dir"
    printf 'string simutask = Solve;\n' > "$case_dir/script/simu.txt"
    rm -rf "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    mkdir -p "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    sed -i -E \
        "s/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5;/; \
         s/^int[[:space:]]+time_integral[[:space:]]*=.*/int time_integral = 1;/; \
         s/^int[[:space:]]+rk_stage[[:space:]]*=.*/int rk_stage = 3;/; \
         s/^Real[[:space:]]+cflst[[:space:]]*=.*/Real cflst = ${cfl};/; \
         s/^Real[[:space:]]+cfled[[:space:]]*=.*/Real cfled = ${cfl};/; \
         s/^int[[:space:]]+ilim[[:space:]]*=.*/int ilim = 0;/" \
        "$case_dir/script/cfd.txt"
    sed -i -E \
        "s/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0;/" \
        "$case_dir/script/vismodel.txt"
    sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
        "$case_dir/script/solver.txt"
    if [ -e "$case_dir/script/grid.txt" ]; then
        sed -i -E \
            "s/^int[[:space:]]+gridObj[[:space:]]*=.*/int gridObj = 3;/" \
            "$case_dir/script/grid.txt"
    fi
}

set_run_controls() {
    local case_dir="$1" max_steps="$2" field_save="$3" start="$4"
    sed -i -E \
        "s/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = ${max_steps};/; \
         s/^int[[:space:]]+nFieldSave[[:space:]]*=.*/int nFieldSave = ${field_save};/; \
         s/^int[[:space:]]+startStrategy[[:space:]]*=.*/int startStrategy = ${start};/" \
        "$case_dir/script/cfd.txt"
}

run_mode() {
    local backend="$1" binary="$2" case_dir="$3" log="$4" breakdown="$5"
    local rc
    set +e
    (
        unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
            ONEFLOW_ENABLE_UNS_HIP_GRADIENT \
            ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
            ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE \
            ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
            ONEFLOW_ENABLE_UNS_HIP_TIME_STEP \
            ONEFLOW_ENABLE_UNS_HIP_INTERFACE || true
        if [ "$backend" = HIP ]; then
            export ONEFLOW_ACCEL_BACKEND=HIP
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
        fi
        export ONEFLOW_FAIL_ON_NONPHYSICAL=1
        export ONEFLOW_STAGE_BREAKDOWN_FILE="$breakdown"
        timeout --foreground --signal=TERM --kill-after=15s \
            "$timeout_seconds" "$binary" 0 "$case_dir"
    ) > "$log" 2>&1
    rc=$?
    set -e
    return "$rc"
}

printf 'backend\tmetric\tvalue\ttolerance\tstatus\n' \
    > "$ONEFLOW_RESTART_RUN_ROOT/restart.tsv"
overall_status=0

for backend in $backends; do
    case "$backend" in
        CPU) binary="$ONEFLOW_RESTART_CPU_BINARY" ;;
        HIP) binary="$ONEFLOW_RESTART_HIP_BINARY" ;;
        *) echo "invalid backend $backend (expected CPU or HIP)" >&2; exit 2 ;;
    esac
    if [ ! -d "$(dirname "$binary")/system" ]; then
        echo "missing $(dirname "$binary")/system" >&2
        exit 2
    fi
    continuous="$ONEFLOW_RESTART_RUN_ROOT/$backend-continuous"
    restarted="$ONEFLOW_RESTART_RUN_ROOT/$backend-restart"

    prepare_case "$continuous"
    set_run_controls "$continuous" "$steps" "$steps" 0
    run_mode "$backend" "$binary" "$continuous" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-continuous.log" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-continuous-%r.tsv"

    prepare_case "$restarted"
    set_run_controls "$restarted" "$split" "$split" 0
    run_mode "$backend" "$binary" "$restarted" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-first-half.log" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-first-half-%r.tsv"
    set_run_controls "$restarted" "$steps" "$steps" 1
    run_mode "$backend" "$binary" "$restarted" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-second-half.log" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-second-half-%r.tsv"

    if [ -r "$rank_verifier" ]; then
        python3 "$rank_verifier" --expected-ranks 1 \
            --label "$backend/continuous" \
            --pattern "$ONEFLOW_RESTART_RUN_ROOT/$backend-continuous-*.tsv" \
            >> "$ONEFLOW_RESTART_RUN_ROOT/rank-utilization.txt"
        python3 "$rank_verifier" --expected-ranks 1 \
            --label "$backend/restart" \
            --pattern "$ONEFLOW_RESTART_RUN_ROOT/$backend-second-half-*.tsv" \
            >> "$ONEFLOW_RESTART_RUN_ROOT/rank-utilization.txt"
    fi

    python3 - "$backend" \
        "$continuous/restart/flow.dat" "$restarted/restart/flow.dat" \
        "$continuous/results/res.dat" "$restarted/results/res.dat" \
        "$continuous/results/aero.dat" "$restarted/results/aero.dat" \
        "$ONEFLOW_RESTART_RUN_ROOT/$backend-continuous.log" \
        "$flow_tolerance" "$residual_tolerance" "$aero_tolerance" \
        "$ONEFLOW_RESTART_RUN_ROOT/restart.tsv" <<'PY'
import pathlib
import re
import struct
import sys

(backend, flow_a, flow_b, res_a, res_b, aero_a, aero_b, log,
 flow_tol, res_tol, aero_tol, out) = sys.argv[1:]
flow_tol = float(flow_tol)
res_tol = float(res_tol)
aero_tol = float(aero_tol)

def grid_extents(text):
    elements = re.search(r"number of elements\s*:\s*(\d+)", text)
    boundary = re.search(r"nBFaces\s*=\s*(\d+)", text)
    if not elements or not boundary:
        raise SystemExit("cannot determine grid extents from %s" % log)
    return int(elements.group(1)) + int(boundary.group(1))

def compare_arrays(path_a, path_b, n_total, tolerance):
    size_a = pathlib.Path(path_a).stat().st_size
    size_b = pathlib.Path(path_b).stat().st_size
    if size_a != size_b:
        return "file_size_mismatch", "size", "0"
    payload = size_a - 20
    if payload <= 0 or payload % (8 * n_total) != 0:
        return "unexpected_layout", str(payload), "0"
    n_arrays = payload // (8 * n_total)
    with open(path_a, "rb") as stream_a, open(path_b, "rb") as stream_b:
        stream_a.read(20)
        stream_b.read(20)
        max_rel = 0.0
        for _ in range(n_arrays):
            a = struct.unpack("<%dd" % n_total, stream_a.read(8 * n_total))
            b = struct.unpack("<%dd" % n_total, stream_b.read(8 * n_total))
            for x, y in zip(a, b):
                if x != y:
                    rel = abs(x - y) / max(1.0, abs(x))
                    if rel > max_rel:
                        max_rel = rel
    return ("PASS" if max_rel <= tolerance else "FAIL"), "%.6e" % max_rel, "%.1e" % tolerance

def compare_last_row(path_a, path_b, tolerance):
    def last(path):
        rows = []
        for line in open(path, errors="replace"):
            fields = line.split()
            if fields and fields[0].lstrip("+-").isdigit():
                rows.append(fields)
        return [float(x) for x in rows[-1][2:]]
    left = last(path_a)
    right = last(path_b)
    if len(left) != len(right):
        return "FAIL", "dimension_mismatch", "%.1e" % tolerance
    max_rel = max(
        abs(x - y) / max(1.0, abs(x)) for x, y in zip(left, right))
    return ("PASS" if max_rel <= tolerance else "FAIL"), "%.6e" % max_rel, "%.1e" % tolerance

log_text = pathlib.Path(log).read_text(errors="replace")
n_total = grid_extents(log_text)
rows = [
    ("flow",) + compare_arrays(flow_a, flow_b, n_total, flow_tol),
    ("residual",) + compare_last_row(res_a, res_b, res_tol),
    ("aero",) + compare_last_row(aero_a, aero_b, aero_tol),
]
with open(out, "a") as stream:
    for metric, status, value, tolerance in rows:
        stream.write("%s\t%s\t%s\t%s\t%s\n" % (
            backend, metric, value, tolerance, status))
PY

    while IFS=$'\t' read -r row_backend metric value tolerance status; do
        if [ "$row_backend" = "$backend" ] && [ "$status" != PASS ]; then
            overall_status=1
        fi
    done < "$ONEFLOW_RESTART_RUN_ROOT/restart.tsv"
done

cat "$ONEFLOW_RESTART_RUN_ROOT/restart.tsv"
if [ "$overall_status" -ne 0 ]; then
    echo "RESTART_EQUIVALENCE_FAIL" >&2
    exit 1
fi
echo "RESTART_EQUIVALENCE_PASS steps=$steps split=$split cfl=$cfl"
