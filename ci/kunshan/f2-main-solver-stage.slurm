#!/bin/bash -l
# Compare main-solver traces for legacy CPU, CPU batch, and HIP batch.  An
# explicit benchmark opt-in can then time the same three paths on this node.
#SBATCH --job-name=oneflow_f2_stage
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

: "${ONEFLOW_F2_STAGE_SOURCE_DIR:?Set ONEFLOW_F2_STAGE_SOURCE_DIR}"
: "${ONEFLOW_F2_STAGE_BINARY:?Set ONEFLOW_F2_STAGE_BINARY}"
: "${ONEFLOW_F2_STAGE_CASE:?Set ONEFLOW_F2_STAGE_CASE (relative to source)}"
: "${ONEFLOW_F2_STAGE_STEPS:?Set ONEFLOW_F2_STAGE_STEPS}"
: "${ONEFLOW_F2_STAGE_TIME_INTEGRAL:?Set ONEFLOW_F2_STAGE_TIME_INTEGRAL}"
: "${ONEFLOW_F2_STAGE_RK_STAGES:?Set ONEFLOW_F2_STAGE_RK_STAGES}"
: "${ONEFLOW_F2_STAGE_CFL:-}"

module purge
module load gcc/16.2.0
module load openmpi/5.0.11
module load cmake/4.4.3
module load compiler/dtk/26.04
module load python/3.8.10
export OMPI_MCA_coll="^hcoll"

source_case="$ONEFLOW_F2_STAGE_SOURCE_DIR/$ONEFLOW_F2_STAGE_CASE"
trace_root="${ONEFLOW_F2_STAGE_RUN_DIR:-$PWD/f2-main-solver-stage-${SLURM_JOB_ID:-manual}}"
trace_timeout="${ONEFLOW_F2_STAGE_TRACE_TIMEOUT_SECONDS:-180}"
rank_utilization_verifier="$ONEFLOW_F2_STAGE_SOURCE_DIR/ci/kunshan/verify-rank-utilization.py"

case "$trace_timeout" in
    ''|*[!0-9]*|0) echo "Invalid ONEFLOW_F2_STAGE_TRACE_TIMEOUT_SECONDS=$trace_timeout" >&2; exit 2 ;;
esac

if [ -e "$trace_root" ]; then
    echo "Refusing to overwrite existing F2 stage run directory" >&2
    exit 2
fi
test -r "$rank_utilization_verifier"
mkdir -p "$trace_root"
{
    echo "gcc=$(gcc --version | head -1)"
    echo "cmake=$(cmake --version | head -1)"
    echo "mpi=$(mpirun --version | head -1)"
    echo "mpirun=$(command -v mpirun)"
    echo "OMPI_MCA_coll=$OMPI_MCA_coll"
    echo "OMPI_MCA_coll_hcoll_enable=${OMPI_MCA_coll_hcoll_enable:-unset}"
    echo "hip_boundary=${ONEFLOW_F2_STAGE_HIP_BOUNDARY:-0}"
    echo "hip_time_step=${ONEFLOW_F2_STAGE_HIP_TIME_STEP:-0}"
} > "$trace_root/environment.txt"
grep -q '16\.2\.0' "$trace_root/environment.txt"
grep -q '4\.4\.3' "$trace_root/environment.txt"
grep -q '5\.0\.11' "$trace_root/environment.txt"
for mode in legacy cpu_batch hip_batch; do
    cp -a "$source_case" "$trace_root/$mode"
    rm -rf "$trace_root/$mode/results" "$trace_root/$mode/restart" \
        "$trace_root/$mode/log"
    mkdir -p "$trace_root/$mode/results" "$trace_root/$mode/restart" \
        "$trace_root/$mode/log"
    sed -i -E "s/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = ${ONEFLOW_F2_STAGE_STEPS};/" \
        "$trace_root/$mode/script/cfd.txt"
    sed -i -E 's/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5; \/\/ Lax-Friedrichs F stage gate/' \
        "$trace_root/$mode/script/cfd.txt"
    sed -i -E "s/^int[[:space:]]+time_integral[[:space:]]*=.*/int time_integral = ${ONEFLOW_F2_STAGE_TIME_INTEGRAL};/" \
        "$trace_root/$mode/script/cfd.txt"
    sed -i -E "s/^int[[:space:]]+rk_stage[[:space:]]*=.*/int rk_stage = ${ONEFLOW_F2_STAGE_RK_STAGES};/" \
        "$trace_root/$mode/script/cfd.txt"
    if [ "${ONEFLOW_F2_STAGE_DISABLE_VISCOUS:-0}" = 1 ]; then
        sed -i -E 's/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0; \/\/ F stage inviscid gate/' \
            "$trace_root/$mode/script/vismodel.txt"
        sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
            "$trace_root/$mode/script/solver.txt"
    fi
    if [ -n "${ONEFLOW_F2_STAGE_CFL:-}" ]; then
        sed -i -E "s/^Real[[:space:]]+cflst[[:space:]]*=.*/Real cflst = ${ONEFLOW_F2_STAGE_CFL};/" "$trace_root/$mode/script/cfd.txt"
        sed -i -E "s/^Real[[:space:]]+cfled[[:space:]]*=.*/Real cfled = ${ONEFLOW_F2_STAGE_CFL};/" "$trace_root/$mode/script/cfd.txt"
    fi
    if [ "${ONEFLOW_F2_STAGE_DISABLE_LIMITER:-0}" = 1 ]; then
        sed -i -E 's/^int[[:space:]]+ilim[[:space:]]*=.*/int ilim = 0; \/\/ F stage limiter gate/' \
            "$trace_root/$mode/script/cfd.txt"
    fi
done

binary_dir=$(dirname "$ONEFLOW_F2_STAGE_BINARY")
if [ ! -d "$binary_dir/system" ]; then
    cp -a "$ONEFLOW_F2_STAGE_SOURCE_DIR/system" "$binary_dir/system"
fi

run_case() {
    local mode=$1
    rm -f "$trace_root/$mode.breakdown-"*.tsv
    rm -f "$trace_root/$mode.stage"
    (
        cd "$trace_root"
        export ONEFLOW_STAGE_BREAKDOWN_FILE="$trace_root/$mode.breakdown-%r.tsv"
        export ONEFLOW_UNS_STAGE_TRACE_FILE="$trace_root/$mode.stage"
        unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
            ONEFLOW_ENABLE_UNS_HIP_GRADIENT ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
            ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
            ONEFLOW_ENABLE_UNS_HIP_TIME_STEP || true
        if [ "$mode" = hip_batch ]; then
            export ONEFLOW_ACCEL_BACKEND=HIP
            for feature in GRADIENT RECONSTRUCTION STATE_UPDATE BOUNDARY TIME_STEP; do
                request_var="ONEFLOW_F2_STAGE_HIP_${feature}"
                value="${!request_var:-}"
                if [ "$value" = 0 ] || [ "$value" = 1 ]; then
                    export "ONEFLOW_ENABLE_UNS_HIP_${feature}=$value"
                fi
            done
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
            if [ "$mode" = cpu_batch ]; then
                export ONEFLOW_ENABLE_UNS_CPU_BATCH=1
            fi
        fi
        timeout --signal=TERM --kill-after=10s "$trace_timeout" \
            mpirun -np 1 "$ONEFLOW_F2_STAGE_BINARY" 0 "$mode"
    ) > "$trace_root/$mode.log" 2>&1
}

run_case legacy
run_case cpu_batch
run_case hip_batch

python3 "$ONEFLOW_F2_STAGE_SOURCE_DIR/ci/kunshan/f2-main-solver-stage-verify-stream.py" \
    "$trace_root/legacy.stage" \
    "$trace_root/cpu_batch.stage" \
    "$trace_root/hip_batch.stage" | tee "$trace_root/verify.txt"

# Performance is deliberately opt-in.  The trace gate above always runs first
# so a timing result cannot be reported without the same-step numerical and
# physical checks passing.  Timing runs do not enable the trace writer because
# trace I/O would dominate the measurement on the 3D m6 case.
if [ "${ONEFLOW_F2_STAGE_BENCHMARK:-0}" = 1 ]; then
    bench_steps="${ONEFLOW_F2_STAGE_BENCH_STEPS:-$ONEFLOW_F2_STAGE_STEPS}"
    bench_repeats="${ONEFLOW_F2_STAGE_BENCH_REPEATS:-3}"
    bench_warmup="${ONEFLOW_F2_STAGE_BENCH_WARMUP:-1}"
    bench_timeout="${ONEFLOW_F2_STAGE_BENCH_TIMEOUT_SECONDS:-180}"
    bench_cpu_ranks="${ONEFLOW_F2_STAGE_BENCH_CPU_RANKS:-1}"
    bench_hip_ranks="${ONEFLOW_F2_STAGE_BENCH_HIP_RANKS:-1}"
    bench_root="$trace_root/benchmark"
    timing_file="$bench_root/timing.tsv"

    case "$bench_steps" in
        ''|*[!0-9]*|0) echo "Invalid ONEFLOW_F2_STAGE_BENCH_STEPS=$bench_steps" >&2; exit 2 ;;
    esac
    case "$bench_repeats" in
        ''|*[!0-9]*|0) echo "Invalid ONEFLOW_F2_STAGE_BENCH_REPEATS=$bench_repeats" >&2; exit 2 ;;
    esac
    case "$bench_warmup" in
        ''|*[!0-9]*) echo "Invalid ONEFLOW_F2_STAGE_BENCH_WARMUP=$bench_warmup" >&2; exit 2 ;;
    esac
    case "$bench_cpu_ranks" in
        ''|*[!0-9]*|0) echo "Invalid ONEFLOW_F2_STAGE_BENCH_CPU_RANKS=$bench_cpu_ranks" >&2; exit 2 ;;
    esac
    case "$bench_hip_ranks" in
        ''|*[!0-9]*|0) echo "Invalid ONEFLOW_F2_STAGE_BENCH_HIP_RANKS=$bench_hip_ranks" >&2; exit 2 ;;
    esac
    if [ "$bench_steps" != "$ONEFLOW_F2_STAGE_STEPS" ]; then
        echo "Benchmark steps must equal trace steps: benchmark=$bench_steps trace=$ONEFLOW_F2_STAGE_STEPS" >&2
        exit 2
    fi

    mkdir -p "$bench_root"
    printf 'mode\trepeat\tranks\telapsed_ms\tsteps\n' > "$timing_file"

    configure_mode() {
        local mode=$1
        unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
            ONEFLOW_ENABLE_UNS_HIP_GRADIENT ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
            ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
            ONEFLOW_ENABLE_UNS_HIP_TIME_STEP ONEFLOW_UNS_STAGE_TRACE_FILE || true
        if [ "$mode" = hip_batch ]; then
            export ONEFLOW_ACCEL_BACKEND=HIP
            for feature in GRADIENT RECONSTRUCTION STATE_UPDATE BOUNDARY TIME_STEP; do
                request_var="ONEFLOW_F2_STAGE_HIP_${feature}"
                value="${!request_var:-}"
                if [ "$value" = 0 ] || [ "$value" = 1 ]; then
                    export "ONEFLOW_ENABLE_UNS_HIP_${feature}=$value"
                fi
            done
        else
            export ONEFLOW_ACCEL_BACKEND=CPU
            if [ "$mode" = cpu_batch ]; then
                export ONEFLOW_ENABLE_UNS_CPU_BATCH=1
            fi
        fi
    }

    reset_benchmark_case() {
        local mode=$1
        rm -rf "$bench_root/$mode/results" "$bench_root/$mode/restart" \
            "$bench_root/$mode/log"
        mkdir -p "$bench_root/$mode/results" "$bench_root/$mode/restart" \
            "$bench_root/$mode/log"
    }

    run_benchmark_case() {
        local mode=$1
        local repeat start_ns end_ns elapsed_ms log_file bench_ranks

        if [ "$mode" = hip_batch ]; then
            bench_ranks="$bench_hip_ranks"
        else
            bench_ranks="$bench_cpu_ranks"
        fi

        for repeat in $(seq 1 "$bench_warmup"); do
            reset_benchmark_case "$mode"
            log_file="$bench_root/$mode/warmup-$repeat.log"
            (
                cd "$bench_root"
                mkdir -p "$bench_root/$mode/warmup-$repeat"
                export ONEFLOW_STAGE_BREAKDOWN_FILE="$bench_root/$mode/warmup-$repeat/breakdown-%r.tsv"
                configure_mode "$mode"
                timeout --signal=TERM --kill-after=10s "$bench_timeout" \
                    mpirun -np "$bench_ranks" "$ONEFLOW_F2_STAGE_BINARY" 0 "$mode"
            ) > "$log_file" 2>&1
            python3 "$rank_utilization_verifier" \
                --expected-ranks "$bench_ranks" \
                --label "$mode/warmup-$repeat" \
                --pattern "$bench_root/$mode/warmup-$repeat/breakdown-*.tsv"
        done

        for repeat in $(seq 1 "$bench_repeats"); do
            reset_benchmark_case "$mode"
            log_file="$bench_root/$mode/repeat-$repeat.log"
            start_ns=$(date +%s%N)
            (
                cd "$bench_root"
                mkdir -p "$bench_root/$mode/repeat-$repeat"
                export ONEFLOW_STAGE_BREAKDOWN_FILE="$bench_root/$mode/repeat-$repeat/breakdown-%r.tsv"
                configure_mode "$mode"
                timeout --signal=TERM --kill-after=10s "$bench_timeout" \
                    mpirun -np "$bench_ranks" "$ONEFLOW_F2_STAGE_BINARY" 0 "$mode"
            ) > "$log_file" 2>&1
            end_ns=$(date +%s%N)
            python3 "$rank_utilization_verifier" \
                --expected-ranks "$bench_ranks" \
                --label "$mode/repeat-$repeat" \
                --pattern "$bench_root/$mode/repeat-$repeat/breakdown-*.tsv"
            elapsed_ms=$(awk -v start="$start_ns" -v end="$end_ns" \
                'BEGIN { printf "%.6f", (end - start) / 1000000.0 }')
            printf '%s\t%s\t%s\t%s\t%s\n' "$mode" "$repeat" \
                "$bench_ranks" "$elapsed_ms" "$bench_steps" >> "$timing_file"
        done
    }

    for mode in legacy cpu_batch hip_batch; do
        rm -rf "$bench_root/$mode"
        cp -a "$trace_root/$mode" "$bench_root/$mode"
        sed -i -E "s/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = ${bench_steps};/" \
            "$bench_root/$mode/script/cfd.txt"
        run_benchmark_case "$mode"
    done

    breakdown_file="$trace_root/breakdown.tsv"
    printf "mode\trepeat\tcategory\telapsed_ms\n" > "$breakdown_file"
    breakdown_count=0
    while IFS= read -r -d "" breakdown_path; do
        breakdown_count=$((breakdown_count + 1))
        repeat_name="$(basename "$(dirname "$breakdown_path")")"
        mode_name="$(basename "$(dirname "$(dirname "$breakdown_path")")")"
        while IFS="$(printf "\t")" read -r category elapsed; do
            printf "%s\t%s\t%s\t%s\n" "$mode_name" "$repeat_name" "$category" "$elapsed" >> "$breakdown_file"
        done < <(tail -n +2 "$breakdown_path")
    done < <(find "$bench_root" -type f -path "*/repeat-*/breakdown-*.tsv" -print0)
    if [ "$breakdown_count" -eq 0 ]; then
        echo "ERROR: no formal-repeat stage breakdown files found" >&2
        exit 1
    fi
    echo "BREAKDOWN_PASS file=$breakdown_file"
    python3 - "$timing_file" <<PY
import collections
import pathlib
import statistics
import sys

path = pathlib.Path(sys.argv[1])
rows = []
for line in path.read_text().splitlines()[1:]:
    mode, repeat, ranks, elapsed_ms, steps = line.split("\t")
    rows.append((mode, int(repeat), int(ranks), float(elapsed_ms), int(steps)))

by_mode = collections.defaultdict(list)
rank_by_mode = {}
for mode, _repeat, ranks, elapsed_ms, _steps in rows:
    by_mode[mode].append(elapsed_ms)
    rank_by_mode.setdefault(mode, ranks)
    if rank_by_mode[mode] != ranks:
        raise SystemExit("inconsistent rank count for mode {}".format(mode))

legacy_mean = statistics.mean(by_mode["legacy"])
print("PERF_BASIS steps={} repeats={} warmup={}".format(
    rows[0][4], len(by_mode["legacy"]), "$bench_warmup"))
for mode in ("legacy", "cpu_batch", "hip_batch"):
    values = by_mode[mode]
    mean = statistics.mean(values)
    median = statistics.median(values)
    speedup = legacy_mean / mean
    print("PERF_SUMMARY mode={} ranks={} mean_ms={:.6f} median_ms={:.6f} "
          "speedup_vs_legacy={:.6f}".format(
              mode, rank_by_mode[mode], mean, median, speedup))
PY
    cp "$timing_file" "$trace_root/timing.tsv"
    echo "PERF_PASS timing_file=$trace_root/timing.tsv"
fi
