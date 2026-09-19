#!/bin/bash -l
# Compare complete main-solver qf1/qf2/invflux traces for legacy and CPU batch.
#SBATCH --job-name=oneflow_e6_cpu_trace
#SBATCH --partition=kshcnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=54G
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

: "${ONEFLOW_E6_TRACE_SOURCE_DIR:?Set ONEFLOW_E6_TRACE_SOURCE_DIR}"
: "${ONEFLOW_E6_TRACE_BINARY:?Set ONEFLOW_E6_TRACE_BINARY}"

module purge
module load compiler/gcc/9.3.0
module load mpi/openmpi/gcc-9.3.0/4.1.5
module load python/3.8.10

source_case="$ONEFLOW_E6_TRACE_SOURCE_DIR/test/m6wingroe_sa"
trace_root="${ONEFLOW_E6_TRACE_RUN_DIR:-$PWD/e6-cpu-trace-${SLURM_JOB_ID:-manual}}"
legacy_case="$trace_root/legacy"
batch_case="$trace_root/batch"

rm -rf "$trace_root"
mkdir -p "$trace_root"
cp -a "$source_case" "$legacy_case"
cp -a "$source_case" "$batch_case"

binary_dir=$(dirname "$ONEFLOW_E6_TRACE_BINARY")
if [ ! -d "$binary_dir/system" ]; then
    mkdir -p "$binary_dir"
    cp -a "$ONEFLOW_E6_TRACE_SOURCE_DIR/system" "$binary_dir/system"
fi

for case_dir in "$legacy_case" "$batch_case"; do
    rm -rf "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    mkdir -p "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    sed -i -E 's/^int[[:space:]]+maxSteps[[:space:]]*=.*/int maxSteps = 1;/'         "$case_dir/script/cfd.txt"
    sed -i -E 's/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme = 5; \/\/ Lax-Friedrichs trace/'         "$case_dir/script/cfd.txt"
done

run_case() {
    local mode=$1
    local case_dir=$2
    (
        cd "$trace_root"
        export ONEFLOW_ACCEL_BACKEND=CPU
        export ONEFLOW_UNS_TRACE_FILE="$trace_root/$mode.trace"
        if [ "$mode" = batch ]; then
            export ONEFLOW_ENABLE_UNS_CPU_BATCH=1
        else
            unset ONEFLOW_ENABLE_UNS_CPU_BATCH || true
        fi
        timeout --signal=TERM --kill-after=10s 180             mpirun -np 1 "$ONEFLOW_E6_TRACE_BINARY" 0 "$mode"
    ) > "$trace_root/$mode.log" 2>&1
}

run_case legacy "$legacy_case"
run_case batch "$batch_case"

python3 "$ONEFLOW_E6_TRACE_SOURCE_DIR/ci/kunshan/e6-cpu-trace-verify.py"     "$trace_root/legacy.trace" "$trace_root/batch.trace"
