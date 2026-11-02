#!/bin/bash
set -euo pipefail

: "${ONEFLOW_MZ_SOURCE_CASE:?Set ONEFLOW_MZ_SOURCE_CASE}"
: "${ONEFLOW_MZ_CPU_BINARY:?Set ONEFLOW_MZ_CPU_BINARY}"
: "${ONEFLOW_MZ_HIP_BINARY:?Set ONEFLOW_MZ_HIP_BINARY}"
: "${ONEFLOW_MZ_RUN_ROOT:?Set ONEFLOW_MZ_RUN_ROOT}"
: "${ONEFLOW_MZ_TRACE_VERIFIER:?Set ONEFLOW_MZ_TRACE_VERIFIER}"
: "${ONEFLOW_MZ_RANK_VERIFIER:?Set ONEFLOW_MZ_RANK_VERIFIER}"

ranks="${ONEFLOW_MZ_RANKS:-8}"
steps="${ONEFLOW_MZ_STEPS:-3}"
cfl="${ONEFLOW_MZ_CFL:-0.01}"
timeout_seconds="${ONEFLOW_MZ_TIMEOUT_SECONDS:-900}"
# MPI placement: `slot` fills nodes in order (fine on one node).  A multi-node
# allocation must pin the per-node rank count, otherwise the scheduler can
# hand out an uneven split (7+1 was observed, which then crashed OpenMPI's
# coll/han communicator setup).
map_by="${ONEFLOW_MZ_MAP_BY:-slot}"
# Stage opt-ins are explicit: the default trace only exercises the fused flux
# and interface path.  ONEFLOW_MZ_HIP_GRADIENT=1 additionally re-opens the
# device gradient, and ONEFLOW_MZ_RESIDENT_CONTRACT=1 verifies that the
# state-owned device q still matches the host mirror before every gradient.
hip_gradient="${ONEFLOW_MZ_HIP_GRADIENT:-}"
hip_reconstruction="${ONEFLOW_MZ_HIP_RECONSTRUCTION:-}"
hip_state_update="${ONEFLOW_MZ_HIP_STATE_UPDATE:-}"
hip_boundary="${ONEFLOW_MZ_HIP_BOUNDARY:-}"
hip_time_step="${ONEFLOW_MZ_HIP_TIME_STEP:-}"
resident_contract="${ONEFLOW_MZ_RESIDENT_CONTRACT:-0}"

for name in ranks steps timeout_seconds; do
    value="${!name}"
    case "$value" in
        ''|*[!0-9]*|0) echo "invalid $name=$value" >&2; exit 2 ;;
    esac
done
if ! awk -v value="$cfl" 'BEGIN {
    valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/;
    exit !(valid && value + 0.0 > 0.0);
}'; then
    echo "invalid cfl=$cfl" >&2
    exit 2
fi
test -d "$ONEFLOW_MZ_SOURCE_CASE"
test -x "$ONEFLOW_MZ_CPU_BINARY"
test -x "$ONEFLOW_MZ_HIP_BINARY"
test -r "$ONEFLOW_MZ_TRACE_VERIFIER"
test -r "$ONEFLOW_MZ_RANK_VERIFIER"
fixture_verifier="${ONEFLOW_MZ_FIXTURE_VERIFIER:-$(dirname "$0")/verify-multizone-fixture.py}"
fixture_merger="${ONEFLOW_MZ_FIXTURE_MERGER:-$(dirname "$0")/merge-multizone-fixture.py}"
fixture_metadata="${ONEFLOW_MZ_FIXTURE_METADATA:-${KUNSHAN_MZ_FIXTURE_METADATA:-}}"
fixture_export_dir="${ONEFLOW_MZ_FIXTURE_EXPORT_DIR:-}"
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
if [ -n "$fixture_export_dir" ]; then
    test -r "$fixture_merger"
    mkdir -p "$fixture_export_dir"
    if compgen -G "$fixture_export_dir/rank-*.json" > /dev/null; then
        echo "fixture export directory already contains rank shards: $fixture_export_dir" >&2
        exit 2
    fi
fi
if [ -e "$ONEFLOW_MZ_RUN_ROOT" ]; then
    echo "refusing to overwrite $ONEFLOW_MZ_RUN_ROOT" >&2
    exit 2
fi
mkdir -p "$ONEFLOW_MZ_RUN_ROOT"
if [ -z "$fixture_metadata" ] && [ -z "$fixture_export_dir" ]; then
    fixture_export_dir="$ONEFLOW_MZ_RUN_ROOT/fixture-shards"
fi
if [ -n "$fixture_export_dir" ]; then
    mkdir -p "$fixture_export_dir"
fi
printf 'fixture_metadata=%s\nfixture_summary=%s\n' \
    "${fixture_metadata:-unset}" "$fixture_summary" \
    > "$ONEFLOW_MZ_RUN_ROOT/fixture-validation.txt"

prepare_case() {
    local case_dir="$1"
    cp -a "$ONEFLOW_MZ_SOURCE_CASE" "$case_dir"
    printf 'string simutask = Solve;\n' > "$case_dir/script/simu.txt"
    printf '4\n1 BCFarfield1\n0 BCOutflow0\n2 BCSymmetryPlane2\n3 BCWall0\n' \
        > "$case_dir/grid/bcRegionMap.txt"
    rm -rf "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    mkdir -p "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    sed -i -E \
        "s/^int[[:space:]]+gridObj[[:space:]]*=.*/int gridObj = 3;/" \
        "$case_dir/script/grid.txt"
    sed -i -E \
        's#^string gridFileName.*#string gridFileName = grid/m6wing-8.ofl;#' \
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
    sed -i -E \
        "s/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0;/; \
         /^[[:space:]]*TurbSolver[[:space:]]*$/d" \
        "$case_dir/script/vismodel.txt"
    sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
        "$case_dir/script/solver.txt"
}

run_mode() {
    local mode="$1"
    local binary="$2"
    local case_dir="$ONEFLOW_MZ_RUN_ROOT/case-$mode"
    local log="$ONEFLOW_MZ_RUN_ROOT/$mode.log"
    local trace="$ONEFLOW_MZ_RUN_ROOT/$mode-%r.stage"
    local breakdown="$ONEFLOW_MZ_RUN_ROOT/$mode-%r.tsv"

    prepare_case "$case_dir"
    (
        cd "$ONEFLOW_MZ_RUN_ROOT"
        unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
            ONEFLOW_ENABLE_UNS_HIP_GRADIENT ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
            ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
            ONEFLOW_ENABLE_UNS_HIP_INTERFACE ONEFLOW_ENABLE_UNS_HIP_TIME_STEP || true
        if [ "$mode" = hip ]; then
            # Only the backend selector is required: the production policy
            # enables the stages that have a passing trace gate (fused flux,
            # resident interface, device gradient).  The legacy stage
            # variables below remain available to force one stage for
            # bisection.
            export ONEFLOW_ACCEL_BACKEND=HIP
            if [ "$hip_gradient" = "1" ]; then
                export ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1
            fi
            if [ "$hip_reconstruction" = "1" ]; then
                export ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION=1
            fi
            if [ "$hip_state_update" = "1" ]; then
                export ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE=1
            fi
            if [ "$hip_boundary" = "1" ]; then
                export ONEFLOW_ENABLE_UNS_HIP_BOUNDARY=1
            fi
            if [ "$hip_time_step" = "1" ]; then
                export ONEFLOW_ENABLE_UNS_HIP_TIME_STEP=1
            fi
            if [ "$resident_contract" = "1" ]; then
                export ONEFLOW_UNS_HIP_RESIDENT_CONTRACT=1
            fi
            export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
            export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
        fi
        export ONEFLOW_STAGE_BREAKDOWN_FILE="$breakdown"
        export ONEFLOW_UNS_STAGE_TRACE_FILE="$trace"
        if [ "$mode" = cpu ] && [ -n "$fixture_export_dir" ]; then
            export ONEFLOW_MZ_FIXTURE_EXPORT_DIR="$fixture_export_dir"
        else
            unset ONEFLOW_MZ_FIXTURE_EXPORT_DIR || true
        fi
        timeout --foreground --signal=TERM --kill-after=15s "$timeout_seconds" \
            mpirun --bind-to core --map-by $map_by --report-bindings -np "$ranks" \
            env ONEFLOW_FAIL_ON_NONPHYSICAL=1 "$binary" 0 "$case_dir"
    ) > "$log" 2>&1
    python3 "$ONEFLOW_MZ_RANK_VERIFIER" \
        --expected-ranks "$ranks" \
        --label "$mode" \
        --pattern "$ONEFLOW_MZ_RUN_ROOT/$mode-*.tsv"
    if [ "$mode" = hip ]; then
        # A CPU fallback must not be reported as HIP evidence: every rank has
        # to show a non-zero fused HIP kernel stage.
        for rank in $(seq 0 $((ranks - 1))); do
            awk -F'\t' '$1 == "HIP kernel" { seen = 1; if ($2 + 0 > 0) ran = 1 }
                END { exit !(seen && ran) }' \
                "$ONEFLOW_MZ_RUN_ROOT/hip-${rank}.tsv" \
                || { echo "HIP kernel stage missing or zero on rank $rank" >&2; exit 1; }
        done
    fi
}

run_mode cpu "$ONEFLOW_MZ_CPU_BINARY"
if [ -n "$fixture_export_dir" ] && [ -z "$fixture_metadata" ]; then
    fixture_metadata="$ONEFLOW_MZ_RUN_ROOT/fixture-metadata.json"
    python3 "$fixture_merger" \
        --shards "$fixture_export_dir" \
        --output "$fixture_metadata"
    fixture_summary="$(
        python3 "$fixture_verifier" \
            --metadata "$fixture_metadata" \
            --expected-ranks "$ranks"
    )"
    printf 'fixture_metadata=%s\nfixture_summary=%s\n' \
        "$fixture_metadata" "$fixture_summary" \
        > "$ONEFLOW_MZ_RUN_ROOT/fixture-validation.txt"
fi
run_mode hip "$ONEFLOW_MZ_HIP_BINARY"

resident_summary="RESIDENT_STATE_CONTRACT_UNCHECKED"
if [ "$resident_contract" = "1" ]; then
    contract_log="$ONEFLOW_MZ_RUN_ROOT/hip.log"
    contract_records="$(grep -c 'RESIDENT_STATE_CONTRACT ' "$contract_log" || true)"
    contract_failures="$(
        grep -c 'RESIDENT_STATE_CONTRACT .*status=FAIL' "$contract_log" || true
    )"
    contract_max_diff="$(
        grep -o 'max_abs_diff=[^ ]*' "$contract_log" \
            | cut -d= -f2 | sort -g | tail -1 || true
    )"
    resident_summary="RESIDENT_STATE_CONTRACT records=${contract_records} \
failures=${contract_failures} max_abs_diff=${contract_max_diff:-unknown}"
    printf '%s\n' "$resident_summary"
    test "$contract_records" -gt 0
    test "$contract_failures" -eq 0
fi
printf 'resident_contract=%s\n' "${resident_summary}" \
    >> "$ONEFLOW_MZ_RUN_ROOT/fixture-validation.txt"

trace_pids=()
for rank in $(seq 0 $((ranks - 1))); do
    python3 "$ONEFLOW_MZ_TRACE_VERIFIER" \
        "$ONEFLOW_MZ_RUN_ROOT/cpu-${rank}.stage" \
        "$ONEFLOW_MZ_RUN_ROOT/cpu-${rank}.stage" \
        "$ONEFLOW_MZ_RUN_ROOT/hip-${rank}.stage" \
        --absolute-tolerance 1e-10 \
        --relative-tolerance 1e-10 \
        > "$ONEFLOW_MZ_RUN_ROOT/trace-rank-${rank}.txt" &
    trace_pids+=("$!")
done
trace_failed=0
for rank in $(seq 0 $((ranks - 1))); do
    if ! wait "${trace_pids[$rank]}"; then
        trace_failed=1
    fi
    grep -q '^STAGE_TRACE_PASS ' "$ONEFLOW_MZ_RUN_ROOT/trace-rank-${rank}.txt" || trace_failed=1
done
test "$trace_failed" -eq 0
printf 'MULTIZONE_TRACE_PASS ranks=%s steps=%s cfl=%s hip_gradient=%s hip_reconstruction=%s hip_state_update=%s hip_boundary=%s hip_time_step=%s resident_contract=%s\n' \
    "$ranks" "$steps" "$cfl" "${hip_gradient:-default}" \
    "${hip_reconstruction:-default}" "${hip_state_update:-default}" \
    "${hip_boundary:-default}" "${hip_time_step:-default}" \
    "$resident_contract" \
    | tee "$ONEFLOW_MZ_RUN_ROOT/result.txt"
