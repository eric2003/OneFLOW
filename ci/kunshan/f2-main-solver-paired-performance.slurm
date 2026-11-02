#!/bin/bash -l
# Paired application-level timing for two validated 3D HIP solver revisions.
# The baseline and candidate HIP runs alternate on one allocated node to reduce
# order and node bias. Accuracy and stability gates must pass before use.
#SBATCH --job-name=oneflow_f2_paired
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

: "${ONEFLOW_F2_PAIRED_SOURCE_DIR:?Set ONEFLOW_F2_PAIRED_SOURCE_DIR}"
: "${ONEFLOW_F2_PAIRED_BASELINE_BINARY:?Set ONEFLOW_F2_PAIRED_BASELINE_BINARY}"
: "${ONEFLOW_F2_PAIRED_CANDIDATE_BINARY:?Set ONEFLOW_F2_PAIRED_CANDIDATE_BINARY}"

if [ "${ONEFLOW_F2_PAIRED_ACCURACY_CONFIRMED:-0}" != 1 ]; then
    echo "Refusing paired timing before the accuracy/stability gates are confirmed" >&2
    exit 2
fi

module purge
module load gcc/16.2.0
module load openmpi/5.0.11
module load cmake/4.4.3
module load compiler/dtk/26.04
module load python/3.8.10

clang_bin="$(command -v clang)"
llvm_root="$(dirname "$(dirname "$clang_bin")")"
dtk_root="$(dirname "$llvm_root")"
mpi_root="$(dirname "$(dirname "$(command -v mpirun)")")"
metis_root="$HOME/.scnet-hpc/oneflow/deps/metis-install"
export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:$dtk_root/dcc/lib64:$dtk_root/dcc/comgr/lib64:$metis_root/lib:$mpi_root/lib:${LD_LIBRARY_PATH:-}"
export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"
export OMPI_MCA_coll="^hcoll"

source_case="$ONEFLOW_F2_PAIRED_SOURCE_DIR/${ONEFLOW_F2_PAIRED_CASE:-test/cfd/m6wingroe_sa}"
run_root="${ONEFLOW_F2_PAIRED_RUN_DIR:-$PWD/f2-main-solver-paired-${SLURM_JOB_ID:-manual}}"
steps="${ONEFLOW_F2_PAIRED_STEPS:-25}"
warmup="${ONEFLOW_F2_PAIRED_WARMUP:-1}"
repeats="${ONEFLOW_F2_PAIRED_REPEATS:-3}"
timeout_seconds="${ONEFLOW_F2_PAIRED_TIMEOUT_SECONDS:-900}"
cpu_ranks="${ONEFLOW_F2_PAIRED_CPU_RANKS:-1}"
hip_ranks="${ONEFLOW_F2_PAIRED_HIP_RANKS:-1}"
cfl="${ONEFLOW_F2_PAIRED_CFL:-0.01}"
time_integral="${ONEFLOW_F2_PAIRED_TIME_INTEGRAL:-1}"
rk_stages="${ONEFLOW_F2_PAIRED_RK_STAGES:-3}"
baseline_label="${ONEFLOW_F2_PAIRED_BASELINE_LABEL:-baseline}"
candidate_label="${ONEFLOW_F2_PAIRED_CANDIDATE_LABEL:-candidate}"
rank_utilization_verifier="$ONEFLOW_F2_PAIRED_SOURCE_DIR/ci/kunshan/verify-rank-utilization.py"

for value_name in steps repeats timeout_seconds cpu_ranks hip_ranks time_integral rk_stages; do
    value="${!value_name}"
    case "$value" in
        ''|*[!0-9]*|0) echo "Invalid $value_name=$value" >&2; exit 2 ;;
    esac
done
case "$warmup" in
    ''|*[!0-9]*) echo "Invalid warmup=$warmup" >&2; exit 2 ;;
esac
if ! awk -v value="$cfl" \
    'BEGIN {
        valid = value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/;
        exit !(valid && value + 0.0 > 0.0);
    }'
then
    echo "Invalid cfl=$cfl" >&2
    exit 2
fi

test -d "$source_case"
test -x "$ONEFLOW_F2_PAIRED_BASELINE_BINARY"
test -x "$ONEFLOW_F2_PAIRED_CANDIDATE_BINARY"
test -r "$rank_utilization_verifier"
if [ -e "$run_root" ]; then
    echo "Refusing to overwrite existing paired run directory" >&2
    exit 2
fi
mkdir -p "$run_root/breakdown"

{
    echo "gcc=$(gcc --version | head -1)"
    echo "cmake=$(cmake --version | head -1)"
    echo "mpi=$(mpirun --version | head -1)"
    echo "mpirun=$(command -v mpirun)"
    echo "OMPI_MCA_coll=$OMPI_MCA_coll"
    echo "hip_arch=gfx906"
    echo "steps=$steps"
    echo "warmup=$warmup"
    echo "repeats=$repeats"
    echo "cpu_ranks=$cpu_ranks"
    echo "hip_ranks=$hip_ranks"
    echo "cfl=$cfl"
    echo "baseline_label=$baseline_label"
    echo "candidate_label=$candidate_label"
} > "$run_root/environment.txt"
grep -q '16\.2\.0' "$run_root/environment.txt"
grep -q '4\.4\.3' "$run_root/environment.txt"
grep -q '5\.0\.11' "$run_root/environment.txt"

prepare_case() {
    local mode=$1
    cp -a "$source_case" "$run_root/$mode"
    rm -rf "$run_root/$mode/results" "$run_root/$mode/restart" \
        "$run_root/$mode/log"
    mkdir -p "$run_root/$mode/results" "$run_root/$mode/restart" \
        "$run_root/$mode/log"
    sed -i -E "s/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = ${steps};/" \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E 's/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5; \/\/ Lax-Friedrichs paired timing/' \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E "s/^int[[:space:]]+time_integral[[:space:]]*=.*/int time_integral = ${time_integral};/" \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E "s/^int[[:space:]]+rk_stage[[:space:]]*=.*/int rk_stage = ${rk_stages};/" \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E "s/^Real[[:space:]]+cflst[[:space:]]*=.*/Real cflst = ${cfl};/" \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E "s/^Real[[:space:]]+cfled[[:space:]]*=.*/Real cfled = ${cfl};/" \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E 's/^int[[:space:]]+ilim[[:space:]]*=.*/int ilim = 0; \/\/ paired timing limiter gate/' \
        "$run_root/$mode/script/cfd.txt"
    sed -i -E 's/^int[[:space:]]+vismodel[[:space:]]*=.*/int vismodel = 0; \/\/ paired timing inviscid gate/' \
        "$run_root/$mode/script/vismodel.txt"
    sed -i -E '/^[[:space:]]*TurbSolver[[:space:]]*$/d' \
        "$run_root/$mode/script/solver.txt"
}

for mode in legacy baseline candidate; do
    prepare_case "$mode"
done

ensure_system_tree() {
    local binary=$1
    local binary_dir
    binary_dir="$(dirname "$binary")"
    if [ ! -d "$binary_dir/system" ]; then
        cp -a "$ONEFLOW_F2_PAIRED_SOURCE_DIR/system" "$binary_dir/system"
    fi
}
ensure_system_tree "$ONEFLOW_F2_PAIRED_BASELINE_BINARY"
ensure_system_tree "$ONEFLOW_F2_PAIRED_CANDIDATE_BINARY"

timing_file="$run_root/timing.tsv"
printf 'phase\tsequence\tmode\trepeat\tranks\telapsed_ms\tsteps\texit_code\n' \
    > "$timing_file"
sequence=0

reset_case() {
    local mode=$1
    rm -rf "$run_root/$mode/results" "$run_root/$mode/restart" \
        "$run_root/$mode/log"
    mkdir -p "$run_root/$mode/results" "$run_root/$mode/restart" \
        "$run_root/$mode/log"
}

configure_mode() {
    local mode=$1
    unset ONEFLOW_ENABLE_UNS_CPU_BATCH ONEFLOW_ENABLE_UNS_HIP_BATCH \
        ONEFLOW_ENABLE_UNS_HIP_GRADIENT ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION \
        ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE ONEFLOW_ENABLE_UNS_HIP_BOUNDARY \
        ONEFLOW_ENABLE_UNS_HIP_TIME_STEP ONEFLOW_UNS_STAGE_TRACE_FILE || true
    if [ "$mode" = legacy ]; then
        export ONEFLOW_ACCEL_BACKEND=CPU
        return
    fi
    export ONEFLOW_ACCEL_BACKEND=HIP
    export ONEFLOW_ENABLE_UNS_HIP_BATCH=1
    export ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1
    export ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION=1
    export ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE=1
    export ONEFLOW_ENABLE_UNS_HIP_BOUNDARY="${ONEFLOW_F2_PAIRED_ENABLE_HIP_BOUNDARY:-0}"
    if [ "$mode" = candidate ]; then
        export ONEFLOW_ENABLE_UNS_HIP_TIME_STEP=1
    fi
}

run_once() {
    local phase=$1
    local mode=$2
    local repeat=$3
    local binary ranks start_ns end_ns elapsed_ms rc log_file

    case "$mode" in
        legacy)
            binary="$ONEFLOW_F2_PAIRED_CANDIDATE_BINARY"
            ranks="$cpu_ranks"
            ;;
        baseline)
            binary="$ONEFLOW_F2_PAIRED_BASELINE_BINARY"
            ranks="$hip_ranks"
            ;;
        candidate)
            binary="$ONEFLOW_F2_PAIRED_CANDIDATE_BINARY"
            ranks="$hip_ranks"
            ;;
        *)
            echo "Unknown timing mode: $mode" >&2
            return 2
            ;;
    esac

    reset_case "$mode"
    sequence=$((sequence + 1))
    log_file="$run_root/${phase}-${sequence}-${mode}-${repeat}.log"
    start_ns="$(date +%s%N)"
    set +e
    (
        cd "$run_root"
        configure_mode "$mode"
        export ONEFLOW_STAGE_BREAKDOWN_FILE="$run_root/breakdown/${phase}-${sequence}-${mode}-${repeat}-%r.tsv"
        timeout --foreground --signal=TERM --kill-after=10s "$timeout_seconds" \
            mpirun -np "$ranks" env ONEFLOW_FAIL_ON_NONPHYSICAL=1 \
            "$binary" 0 "$mode"
    ) > "$log_file" 2>&1
    rc=$?
    set -e
    end_ns="$(date +%s%N)"
    if [ "$rc" -eq 0 ] && ! python3 "$rank_utilization_verifier" \
        --expected-ranks "$ranks" \
        --label "$phase/$sequence/$mode/$repeat" \
        --pattern "$run_root/breakdown/${phase}-${sequence}-${mode}-${repeat}-*.tsv"
    then
        rc=3
    fi
    elapsed_ms="$(awk -v start="$start_ns" -v end="$end_ns" \
        'BEGIN { printf "%.6f", (end - start) / 1000000.0 }')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$phase" "$sequence" "$mode" "$repeat" "$ranks" "$elapsed_ms" \
        "$steps" "$rc" | tee -a "$timing_file"
    return "$rc"
}

for repeat in $(seq 1 "$warmup"); do
    run_once warmup legacy "$repeat"
    if [ $((repeat % 2)) -eq 1 ]; then
        run_once warmup baseline "$repeat"
        run_once warmup candidate "$repeat"
    else
        run_once warmup candidate "$repeat"
        run_once warmup baseline "$repeat"
    fi
done

for repeat in $(seq 1 "$repeats"); do
    run_once formal legacy "$repeat"
done
for repeat in $(seq 1 "$repeats"); do
    if [ $((repeat % 2)) -eq 1 ]; then
        run_once formal baseline "$repeat"
        run_once formal candidate "$repeat"
    else
        run_once formal candidate "$repeat"
        run_once formal baseline "$repeat"
    fi
done

python3 - "$timing_file" "$baseline_label" "$candidate_label" <<'PY'
import collections
import pathlib
import statistics
import sys

path = pathlib.Path(sys.argv[1])
baseline_label = sys.argv[2]
candidate_label = sys.argv[3]
values = collections.defaultdict(list)
for line in path.read_text().splitlines()[1:]:
    phase, _sequence, mode, _repeat, _ranks, elapsed_ms, _steps, exit_code = (
        line.split("\t")
    )
    if phase == "formal":
        if exit_code != "0":
            raise SystemExit("formal timing contains a nonzero workload exit")
        values[mode].append(float(elapsed_ms))

legacy = statistics.mean(values["legacy"])
baseline = statistics.mean(values["baseline"])
candidate = statistics.mean(values["candidate"])
print(
    "PAIRED_BASIS repeats={} legacy_mean_ms={:.6f} "
    "{}_mean_ms={:.6f} {}_mean_ms={:.6f}".format(
        len(values["candidate"]),
        legacy,
        baseline_label,
        baseline,
        candidate_label,
        candidate,
    )
)
print(
    "PAIRED_SUMMARY legacy_over_baseline={:.6f} "
    "legacy_over_candidate={:.6f} baseline_over_candidate={:.6f}".format(
        legacy / baseline,
        legacy / candidate,
        baseline / candidate,
    )
)
PY

printf 'result=0\n' > "$run_root/result.txt"
