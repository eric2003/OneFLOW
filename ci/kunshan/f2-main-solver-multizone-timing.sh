#!/bin/bash
set -euo pipefail

: "${ONEFLOW_MZ_SOURCE_CASE:?Set ONEFLOW_MZ_SOURCE_CASE}"
: "${ONEFLOW_MZ_BINARY:?Set ONEFLOW_MZ_BINARY}"
: "${ONEFLOW_MZ_RUN_ROOT:?Set ONEFLOW_MZ_RUN_ROOT}"
: "${ONEFLOW_MZ_MODE:?Set ONEFLOW_MZ_MODE to cpu or hip}"

mode="$ONEFLOW_MZ_MODE"
ranks="${ONEFLOW_MZ_RANKS:-8}"
steps="${ONEFLOW_MZ_STEPS:-10}"
warmup="${ONEFLOW_MZ_WARMUP:-1}"
repeats="${ONEFLOW_MZ_REPEATS:-5}"
cfl="${ONEFLOW_MZ_CFL:-0.01}"
timeout_seconds="${ONEFLOW_MZ_TIMEOUT_SECONDS:-600}"
# MPI placement: `slot` fills nodes in order (fine on one node).  A multi-node
# allocation must pin the per-node rank count, otherwise the scheduler can
# hand out an uneven split (7+1 was observed, which crashed OpenMPI's
# coll/han communicator setup during case setup).
map_by="${ONEFLOW_MZ_MAP_BY:-slot}"
# The partitioned fixture name is a property of the case, not of the runner.
# Keep the validated 8-zone default and allow other prepared fixtures (for
# example the 4-zone mesh) to reuse the same preparation and timing contract.
grid_file="${ONEFLOW_MZ_GRID_FILE:-grid/m6wing-8.ofl}"
# Optional output cadence.  Leave unset to keep the case defaults; set all three
# above the step count to measure the solver path without intermediate restart
# or aerodynamic files (the final postprocess still runs once).
field_save="${ONEFLOW_MZ_FIELD_SAVE:-}"
res_save="${ONEFLOW_MZ_RES_SAVE:-}"
force_save="${ONEFLOW_MZ_FORCE_SAVE:-}"

case "$mode" in
    cpu|hip) ;;
    *) echo "invalid ONEFLOW_MZ_MODE=$mode" >&2; exit 2 ;;
esac
for name in ranks steps warmup repeats timeout_seconds; do
    value="${!name}"
    case "$value" in
        ''|*[!0-9]*) echo "invalid $name=$value" >&2; exit 2 ;;
    esac
    if [ "$name" != warmup ] && [ "$value" -eq 0 ]; then
        echo "invalid $name=$value" >&2
        exit 2
    fi
done
if ! awk -v value="$cfl" 'BEGIN {
    valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/;
    exit !(valid && value + 0.0 > 0.0);
}'; then
    echo "invalid cfl=$cfl" >&2
    exit 2
fi

rank_verifier="${ONEFLOW_MZ_RANK_VERIFIER:-$(dirname "$ONEFLOW_MZ_BINARY")/verify-rank-utilization.py}"
fixture_verifier="${ONEFLOW_MZ_FIXTURE_VERIFIER:-$(dirname "$0")/verify-multizone-fixture.py}"
test -d "$ONEFLOW_MZ_SOURCE_CASE"
test -x "$ONEFLOW_MZ_BINARY"
test -r "$rank_verifier"
fixture_metadata="${ONEFLOW_MZ_FIXTURE_METADATA:-${KUNSHAN_MZ_FIXTURE_METADATA:-}}"
fixture_summary="MULTIZONE_FIXTURE_UNCHECKED metadata=unset"
if [ -n "$fixture_metadata" ]; then
    test -r "$fixture_metadata"
    test -r "$fixture_verifier"
    fixture_summary="$(
        python3 "$fixture_verifier" \
            --metadata "$fixture_metadata" \
            --expected-ranks "$ranks"
    )"
    printf '%s\n' "$fixture_summary"
fi
if [ -e "$ONEFLOW_MZ_RUN_ROOT" ]; then
    echo "refusing to overwrite $ONEFLOW_MZ_RUN_ROOT" >&2
    exit 2
fi
mkdir -p "$ONEFLOW_MZ_RUN_ROOT"
printf 'fixture_metadata=%s\nfixture_summary=%s\n' \
    "${fixture_metadata:-unset}" "$fixture_summary" \
    > "$ONEFLOW_MZ_RUN_ROOT/fixture-validation.txt"

{
    echo "mode=$mode"
    echo "ranks=$ranks"
    echo "cpu_cores=$ranks"
    echo "steps=$steps"
    echo "warmup=$warmup"
    echo "repeats=$repeats"
    echo "cfl=$cfl"
    echo "dcu_count=$([ "$mode" = hip ] && echo 1 || echo 0)"
    echo "ranks_per_dcu=$([ "$mode" = hip ] && echo "$ranks" || echo 0)"
    echo "fixture_metadata=${fixture_metadata:-unset}"
    echo "fixture_summary=$fixture_summary"
    echo "nodes=${SLURM_JOB_NUM_NODES:-1}"
    echo "hostlist=${SLURM_JOB_NODELIST:-unknown}"
    echo "mpirun=$(command -v mpirun)"
    echo "mpi=$(mpirun --version | head -1)"
    echo "map_by=$map_by"
    echo "grid_file=$grid_file"
    echo "field_save=${field_save:-case-default} res_save=${res_save:-case-default} force_save=${force_save:-case-default}"
    echo "OMPI_MCA_coll=${OMPI_MCA_coll:-unset}"
    echo "OMPI_MCA_coll_hcoll_enable=${OMPI_MCA_coll_hcoll_enable:-unset}"
    echo "SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE:-unknown}"
    echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unset}"
} > "$ONEFLOW_MZ_RUN_ROOT/environment.txt"

timing="$ONEFLOW_MZ_RUN_ROOT/timing.tsv"
printf 'phase\trepeat\tmode\tranks\telapsed_ms\texit_code\n' > "$timing"

prepare_case() {
    local case_dir="$1"
    cp -a "$ONEFLOW_MZ_SOURCE_CASE" "$case_dir"
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
         s/^int[[:space:]]+ilim[[:space:]]*=.*/int ilim = 0;/" \
        "$case_dir/script/cfd.txt"
    if [ -n "$field_save" ]; then
        sed -i -E "s/^int[[:space:]]+nFieldSave[[:space:]]*=.*/int nFieldSave = ${field_save};/" \
            "$case_dir/script/cfd.txt"
    fi
    if [ -n "$res_save" ]; then
        sed -i -E "s/^int[[:space:]]+nResSave[[:space:]]*=.*/int nResSave = ${res_save};/" \
            "$case_dir/script/cfd.txt"
    fi
    if [ -n "$force_save" ]; then
        sed -i -E "s/^int[[:space:]]+nForceSave[[:space:]]*=.*/int nForceSave = ${force_save};/" \
            "$case_dir/script/cfd.txt"
    fi
    sed -i -E \
        "s/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0;/; \
         /^[[:space:]]*TurbSolver[[:space:]]*$/d" \
        "$case_dir/script/vismodel.txt"
    sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
        "$case_dir/script/solver.txt"
}

run_once() {
    local phase="$1"
    local repeat="$2"
    local case_dir="$ONEFLOW_MZ_RUN_ROOT/case-${phase}-${repeat}"
    local log="$ONEFLOW_MZ_RUN_ROOT/${phase}-${repeat}.log"
    local breakdown="$ONEFLOW_MZ_RUN_ROOT/${phase}-${repeat}-breakdown-%r.tsv"
    local start_ns end_ns elapsed rc

    prepare_case "$case_dir"
    start_ns="$(date +%s%N)"
    set +e
    (
        cd "$ONEFLOW_MZ_RUN_ROOT"
        if [ "$mode" = hip ]; then
            # Backend-only: the production policy selects the validated stages.
            export ONEFLOW_ACCEL_BACKEND=HIP
            export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
            export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
            unset ONEFLOW_ENABLE_UNS_HIP_BATCH ONEFLOW_ENABLE_UNS_HIP_GRADIENT \
                ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE \
                ONEFLOW_ENABLE_UNS_HIP_BOUNDARY ONEFLOW_ENABLE_UNS_HIP_INTERFACE || true
        fi
        export ONEFLOW_STAGE_BREAKDOWN_FILE="$breakdown"
        timeout --foreground --signal=TERM --kill-after=15s "$timeout_seconds" \
            mpirun --bind-to core --map-by $map_by --report-bindings -np "$ranks" \
            env ONEFLOW_FAIL_ON_NONPHYSICAL=1 "$ONEFLOW_MZ_BINARY" 0 "$case_dir"
    ) > "$log" 2>&1
    rc=$?
    set -e
    end_ns="$(date +%s%N)"
    if [ "$rc" -eq 0 ]; then
        python3 "$rank_verifier" \
            --expected-ranks "$ranks" \
            --label "$phase/$repeat" \
            --pattern "$ONEFLOW_MZ_RUN_ROOT/${phase}-${repeat}-breakdown-*.tsv" \
            >> "$ONEFLOW_MZ_RUN_ROOT/rank-utilization.txt"
    fi
    elapsed="$(awk -v start="$start_ns" -v end="$end_ns" \
        'BEGIN { printf "%.6f", (end - start) / 1000000.0 }')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$phase" "$repeat" "$mode" "$ranks" "$elapsed" "$rc" | tee -a "$timing"
    if [ "$rc" -ne 0 ]; then
        return "$rc"
    fi
}

for repeat in $(seq 1 "$warmup"); do
    run_once warmup "$repeat"
done
for repeat in $(seq 1 "$repeats"); do
    run_once formal "$repeat"
done

python3 - "$timing" "$ONEFLOW_MZ_RUN_ROOT/summary.txt" <<'PY'
import pathlib
import statistics
import sys

timing = pathlib.Path(sys.argv[1])
summary = pathlib.Path(sys.argv[2])
rows = []
for line in timing.read_text().splitlines()[1:]:
    phase, repeat, mode, ranks, elapsed, rc = line.split("\t")
    if phase == "formal" and rc == "0":
        rows.append(float(elapsed))
if not rows:
    raise SystemExit("no successful formal repeats")
summary.write_text(
    "formal_count=%d\nformal_mean_ms=%.6f\nformal_min_ms=%.6f\n"
    "formal_max_ms=%.6f\n"
    % (len(rows), statistics.mean(rows), min(rows), max(rows))
)
print(summary.read_text(), end="")
PY
