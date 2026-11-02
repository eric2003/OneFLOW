#!/bin/bash
# Cross-backend residual/aero gate for the production 3D solver.
#
# Runs the same partitioned case with the legacy CPU backend and the HIP
# backend for N steps without a trace, then compares the full residual history,
# the final aerodynamic row and (when textual) the final flow file.  This is
# the accuracy gate for host/interface work that changes where a value is
# packed or exchanged: the 3-step trace oracle checks only three steps, while a
# 200-step run catches a drift that accumulates.
#
# Environment:
#   ONEFLOW_XB_SOURCE_CASE   prepared multi-zone case directory (required)
#   ONEFLOW_XB_BINARY        root binary exporting both backends (required)
#   ONEFLOW_XB_RUN_ROOT      append-only output directory (required)
#   ONEFLOW_XB_STEPS         step count (default 200)
#   ONEFLOW_XB_CPU_RANKS     CPU ranks (default 8)
#   ONEFLOW_XB_HIP_RANKS     HIP ranks (default 8, i.e. the 8:1 slot)
#   ONEFLOW_XB_CFL           fixed CFL for both ramp ends (default 0.01)
#   ONEFLOW_XB_GRID_FILE     grid file inside the case (default grid/m6wing-8.ofl)
#   ONEFLOW_XB_TIMEOUT_SECONDS  per-run timeout (default 1800)
#   ONEFLOW_XB_RES_REL_TOL   residual tolerance (default 1e-12)
#   ONEFLOW_XB_AERO_REL_TOL  aero tolerance (default 1e-10)
#   ONEFLOW_XB_FLOW_REL_TOL  flow tolerance (default 1e-12)
#   ONEFLOW_XB_MAP_BY        mpirun map-by (default slot)
#   ONEFLOW_XB_HIP_SKIP_POST_INTERFACE  pass ONEFLOW_SKIP_POST_INTERFACE to the
#                            HIP run only (default unset: legacy behaviour)
#   ONEFLOW_XB_HIP_SKIP_HOST_INTERFACE  pass ONEFLOW_HIP_SKIP_HOST_INTERFACE to
#                            the HIP run only (default unset: legacy behaviour)
#
# The runner intentionally contains no cluster paths, hostnames or job
# metadata. Scope: one local zone per rank, MULTI_STAGE + ilim=0 + vismodel=0,
# which is the validated inviscid device path.
set -euo pipefail

: "${ONEFLOW_XB_SOURCE_CASE:?Set ONEFLOW_XB_SOURCE_CASE}"
: "${ONEFLOW_XB_BINARY:?Set ONEFLOW_XB_BINARY}"
: "${ONEFLOW_XB_RUN_ROOT:?Set ONEFLOW_XB_RUN_ROOT}"

steps="${ONEFLOW_XB_STEPS:-200}"
cpu_ranks="${ONEFLOW_XB_CPU_RANKS:-8}"
hip_ranks="${ONEFLOW_XB_HIP_RANKS:-8}"
cfl="${ONEFLOW_XB_CFL:-0.01}"
grid_file="${ONEFLOW_XB_GRID_FILE:-grid/m6wing-8.ofl}"
timeout_seconds="${ONEFLOW_XB_TIMEOUT_SECONDS:-1800}"
res_tolerance="${ONEFLOW_XB_RES_REL_TOL:-1e-12}"
aero_tolerance="${ONEFLOW_XB_AERO_REL_TOL:-1e-10}"
flow_tolerance="${ONEFLOW_XB_FLOW_REL_TOL:-1e-12}"
map_by="${ONEFLOW_XB_MAP_BY:-slot}"
rank_verifier="${ONEFLOW_XB_RANK_VERIFIER:-$(dirname "$0")/verify-rank-utilization.py}"

for name in steps cpu_ranks hip_ranks timeout_seconds; do
    value="${!name}"
    case "$value" in
        ''|*[!0-9]*|0) echo "invalid $name=$value" >&2; exit 2 ;;
    esac
done
if ! awk -v value="$cfl" 'BEGIN {
    valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/;
    exit !(valid && value + 0.0 > 0.0);
}'; then
    echo "invalid ONEFLOW_XB_CFL=$cfl" >&2
    exit 2
fi

test -d "$ONEFLOW_XB_SOURCE_CASE"
test -x "$ONEFLOW_XB_BINARY"
test -r "$rank_verifier"
if [ -e "$ONEFLOW_XB_RUN_ROOT" ]; then
    echo "refusing to overwrite $ONEFLOW_XB_RUN_ROOT" >&2
    exit 2
fi
mkdir -p "$ONEFLOW_XB_RUN_ROOT"

{
    echo "steps=$steps"
    echo "cpu_ranks=$cpu_ranks"
    echo "hip_ranks=$hip_ranks"
    echo "cfl=$cfl"
    echo "grid_file=$grid_file"
    echo "map_by=$map_by"
    echo "nodes=${SLURM_JOB_NUM_NODES:-1}"
    echo "hostlist=${SLURM_JOB_NODELIST:-unknown}"
    echo "mpi=$(mpirun --version | head -1)"
    echo "OMPI_MCA_coll=${OMPI_MCA_coll:-unset}"
} > "$ONEFLOW_XB_RUN_ROOT/environment.txt"

prepare_case() {
    local case_dir="$1"
    cp -a "$ONEFLOW_XB_SOURCE_CASE" "$case_dir"
    printf 'string simutask = Solve;\n' > "$case_dir/script/simu.txt"
    printf '4\n1 BCFarfield1\n0 BCOutflow0\n2 BCSymmetryPlane2\n3 BCWall0\n' \
        > "$case_dir/grid/bcRegionMap.txt"
    rm -rf "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    mkdir -p "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    if [ -e "$case_dir/script/grid.txt" ]; then
        sed -i -E \
            "s/^int[[:space:]]+gridObj[[:space:]]*=.*/int gridObj = 3;/" \
            "$case_dir/script/grid.txt"
    fi
    sed -i -E \
        "s#^string gridFileName.*#string gridFileName = ${grid_file};#" \
        "$case_dir/script/cfd.txt"
    sed -i -E \
        "s/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = ${steps};/; \
         s/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5;/; \
         s/^int[[:space:]]+time_integral[[:space:]]*=.*/int time_integral = 1;/; \
         s/^int[[:space:]]+rk_stage[[:space:]]*=.*/int rk_stage = 3;/; \
         s/^Real[[:space:]]+cflst[[:space:]]*=.*/Real cflst = ${cfl};/; \
         s/^Real[[:space:]]+cfled[[:space:]]*=.*/Real cfled = ${cfl};/; \
         s/^int[[:space:]]+ilim[[:space:]]*=.*/int ilim = 0;/; \
         s/^int[[:space:]]+nFieldSave[[:space:]]*=.*/int nFieldSave = ${steps};/; \
         s/^int[[:space:]]+nResSave[[:space:]]*=.*/int nResSave = 1;/; \
         s/^int[[:space:]]+nForceSave[[:space:]]*=.*/int nForceSave = 1;/" \
        "$case_dir/script/cfd.txt"
    sed -i -E \
        "s/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0;/; \
         /^[[:space:]]*TurbSolver[[:space:]]*$/d" \
        "$case_dir/script/vismodel.txt"
    sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
        "$case_dir/script/solver.txt"
}

run_backend() {
    local mode="$1" ranks="$2" tag="$3"
    local case_dir="$ONEFLOW_XB_RUN_ROOT/case-$tag"
    local log="$ONEFLOW_XB_RUN_ROOT/$tag.log"
    local rc

    prepare_case "$case_dir"
    set +e
    (
        cd "$ONEFLOW_XB_RUN_ROOT"
        unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
            ONEFLOW_ENABLE_UNS_HIP_GRADIENT ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
            ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
            ONEFLOW_ENABLE_UNS_HIP_TIME_STEP ONEFLOW_ENABLE_UNS_HIP_INTERFACE || true
        # High-precision residual rows: the normal writer rounds too much for a
        # 1e-12 comparison.
        export ONEFLOW_RESIDUAL_TEST_OUTPUT=1
        if [ "$mode" = HIP ]; then
            export ONEFLOW_ACCEL_BACKEND=HIP
            # The CPU run always keeps the legacy host exchange so the
            # comparison is against the unmodified oracle.
            if [ -n "${ONEFLOW_XB_HIP_SKIP_POST_INTERFACE:-}" ]; then
                export ONEFLOW_SKIP_POST_INTERFACE="$ONEFLOW_XB_HIP_SKIP_POST_INTERFACE"
            fi
            if [ -n "${ONEFLOW_XB_HIP_SKIP_HOST_INTERFACE:-}" ]; then
                export ONEFLOW_HIP_SKIP_HOST_INTERFACE="$ONEFLOW_XB_HIP_SKIP_HOST_INTERFACE"
            fi
            if [ -n "${ONEFLOW_XB_HIP_CACHE_HANDSHAKE:-}" ]; then
                export ONEFLOW_HIP_CACHE_INTERFACE_HANDSHAKE="$ONEFLOW_XB_HIP_CACHE_HANDSHAKE"
            fi
            if [ -n "${ONEFLOW_XB_HIP_SKIP_STATE_HOST_SYNC:-}" ]; then
                export ONEFLOW_HIP_SKIP_STATE_HOST_SYNC="$ONEFLOW_XB_HIP_SKIP_STATE_HOST_SYNC"
            fi
            if [ -n "${ONEFLOW_XB_HIP_SKIP_UNUSED_INTERFACE:-}" ]; then
                export ONEFLOW_HIP_SKIP_UNUSED_INTERFACE="$ONEFLOW_XB_HIP_SKIP_UNUSED_INTERFACE"
            fi
            export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
            export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
        fi
        timeout --foreground --signal=TERM --kill-after=15s "$timeout_seconds" \
            mpirun --bind-to core --map-by "$map_by" --report-bindings \
            -np "$ranks" env ONEFLOW_FAIL_ON_NONPHYSICAL=1 \
            "$ONEFLOW_XB_BINARY" 0 "$case_dir"
    ) > "$log" 2>&1
    rc=$?
    set -e

    printf '%s\t%s\t%s\n' "$tag" "$ranks" "$rc" >> "$ONEFLOW_XB_RUN_ROOT/runs.tsv"
    if [ "$rc" -ne 0 ]; then
        echo "CROSS_BACKEND_FAIL backend=$mode rc=$rc" >&2
        return "$rc"
    fi
    test -s "$case_dir/results/res.dat"
    test -s "$case_dir/results/aero.dat"
}

printf 'backend\tranks\texit_code\n' > "$ONEFLOW_XB_RUN_ROOT/runs.tsv"
run_backend CPU "$cpu_ranks" cpu
run_backend HIP "$hip_ranks" hip

python3 - \
    "$ONEFLOW_XB_RUN_ROOT/case-cpu" \
    "$ONEFLOW_XB_RUN_ROOT/case-hip" \
    "$res_tolerance" "$aero_tolerance" "$flow_tolerance" \
    "$ONEFLOW_XB_RUN_ROOT/cross-backend.tsv" <<'PY'
import hashlib
import math
import pathlib
import sys

cpu_dir, hip_dir, res_tol, aero_tol, flow_tol, out_path = sys.argv[1:]
res_tol, aero_tol, flow_tol = float(res_tol), float(aero_tol), float(flow_tol)


def rel(a, b):
    return abs(a - b) / max(1.0, abs(a), abs(b))


def data_rows(path):
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if fields and fields[0].lstrip("+-").isdigit():
            rows.append(fields)
    return rows


def compare_rows(rows_a, rows_b, tol):
    """Compare every numeric field of two equally sized row sets."""
    worst = 0.0
    worst_at = (0, 0)
    for i, (ra, rb) in enumerate(zip(rows_a, rows_b)):
        if len(ra) != len(rb):
            return math.inf, (i, -1)
        for j, (fa, fb) in enumerate(zip(ra, rb)):
            if j == 0:
                continue  # step index
            try:
                va, vb = float(fa), float(fb)
            except ValueError:
                if fa != fb:
                    return math.inf, (i, j)
                continue
            d = rel(va, vb)
            if d > worst:
                worst, worst_at = d, (i, j)
    return worst, worst_at


def flow_compare(path_a, path_b, tol):
    try:
        text_a = path_a.read_text()
        text_b = path_b.read_text()
    except UnicodeDecodeError:
        same = path_a.read_bytes() == path_b.read_bytes()
        return ("binary_byte_equal" if same else "binary_byte_diff"), None
    a = text_a.split()
    b = text_b.split()
    if len(a) != len(b):
        return "dimension_mismatch", None
    worst = 0.0
    for fa, fb in zip(a, b):
        try:
            va, vb = float(fa), float(fb)
        except ValueError:
            if fa != fb:
                return "token_mismatch", None
            continue
        worst = max(worst, rel(va, vb))
    return ("pass" if worst <= tol else "fail"), worst


cpu = pathlib.Path(cpu_dir)
hip = pathlib.Path(hip_dir)

lines = ["check\trows\tmax_rel_diff\trow\tfield\tstatus"]
ok = True

res_a, res_b = data_rows(cpu / "results/res.dat"), data_rows(hip / "results/res.dat")
if len(res_a) != len(res_b) or not res_a:
    lines.append("residual\t%d/%d\tnan\t-\t-\tfail" % (len(res_a), len(res_b)))
    ok = False
else:
    worst, (row, field) = compare_rows(res_a, res_b, res_tol)
    status = "pass" if worst <= res_tol else "fail"
    ok = ok and status == "pass"
    lines.append("residual\t%d\t%.6e\t%d\t%d\t%s"
                 % (len(res_a), worst, row, field, status))

aero_a, aero_b = data_rows(cpu / "results/aero.dat"), data_rows(hip / "results/aero.dat")
if not aero_a or not aero_b:
    lines.append("aero\t%d/%d\tnan\t-\t-\tfail" % (len(aero_a), len(aero_b)))
    ok = False
else:
    worst, (_, field) = compare_rows(aero_a[-1:], aero_b[-1:], aero_tol)
    status = "pass" if worst <= aero_tol else "fail"
    ok = ok and status == "pass"
    lines.append("aero\t%d\t%.6e\tlast\t%d\t%s"
                 % (len(aero_a), worst, field, status))

flow_status, flow_worst = flow_compare(cpu / "restart/flow.dat", hip / "restart/flow.dat",
                                       flow_tol)
if flow_status in ("pass", "binary_byte_equal"):
    pass
elif flow_status == "fail":
    ok = False
elif flow_status in ("dimension_mismatch", "token_mismatch"):
    ok = False
lines.append("flow\t-\t%s\t-\t-\t%s"
             % ("nan" if flow_worst is None else "%.6e" % flow_worst, flow_status))

pathlib.Path(out_path).write_text("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if ok else 1)
PY

echo "CROSS_BACKEND_RESIDUAL_PASS"
