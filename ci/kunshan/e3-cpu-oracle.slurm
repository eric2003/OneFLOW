#!/bin/bash -l
# Compare the legacy and CPU-batch 3D inviscid flux paths on one identical case.
# The source tree and already-built solver are supplied by environment variables.
#SBATCH --job-name=oneflow_e3_cpu_oracle
#SBATCH --partition=kshcnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=54G
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

: "${ONEFLOW_E3_ORACLE_SOURCE_DIR:?Set ONEFLOW_E3_ORACLE_SOURCE_DIR}"
: "${ONEFLOW_E3_ORACLE_BINARY:?Set ONEFLOW_E3_ORACLE_BINARY}"

module purge
module load compiler/gcc/9.3.0
module load mpi/openmpi/gcc-9.3.0/4.1.5
module load python/3.8.10

source_case="$ONEFLOW_E3_ORACLE_SOURCE_DIR/test/m6wingroe_sa"
oracle_root="${ONEFLOW_E3_ORACLE_RUN_DIR:-$PWD/e3-cpu-oracle-${SLURM_JOB_ID:-manual}}"
legacy_case="$oracle_root/legacy/m6wingroe_sa"
batch_case="$oracle_root/batch/m6wingroe_sa"

rm -rf "$oracle_root"
mkdir -p "$oracle_root/legacy" "$oracle_root/batch"
cp -a "$source_case" "$legacy_case"
cp -a "$source_case" "$batch_case"

for case_dir in "$legacy_case" "$batch_case"; do
    rm -rf "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    mkdir -p "$case_dir/results" "$case_dir/restart" "$case_dir/log"
    sed -i -E 's/^int[[:space:]]+ischeme[[:space:]]*=.*/int ischeme   = 5; \/\/ Lax-Friedrichs oracle/' \
        "$case_dir/script/cfd.txt"
done

run_case() {
    local mode=$1
    local case_dir=$2
    (
        cd "$(dirname "$case_dir")"
        export ONEFLOW_ACCEL_BACKEND=CPU
        if [ "$mode" = batch ]; then
            export ONEFLOW_ENABLE_UNS_CPU_BATCH=1
        else
            unset ONEFLOW_ENABLE_UNS_CPU_BATCH || true
        fi
        timeout --signal=TERM --kill-after=10s 180 \
            mpirun -np 1 "$ONEFLOW_E3_ORACLE_BINARY" 0 "$(basename "$case_dir")"
    ) > "$oracle_root/$mode.log" 2>&1
}

run_case legacy "$legacy_case"
run_case batch "$batch_case"

python3 - "$legacy_case" "$batch_case" <<'PY'
import math
import pathlib
import sys

legacy = pathlib.Path(sys.argv[1])
batch = pathlib.Path(sys.argv[2])
files = [
    "results/aero.dat",
    "results/wallaero.dat",
    "results/res.dat",
    "results/turbres.dat",
]
absolute_limit = 2.0e-12
relative_limit = 1.0e-10
overall_abs = 0.0
overall_rel = 0.0

def numbers(path):
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        for token in line.split():
            try:
                yield float(token)
            except ValueError:
                continue

for rel in files:
    left = list(numbers(legacy / rel))
    right = list(numbers(batch / rel))
    if len(left) != len(right):
        raise SystemExit(f"ORACLE_FAIL file={rel} count={len(left)} vs {len(right)}")
    file_abs = 0.0
    file_rel = 0.0
    for a, b in zip(left, right):
        delta = abs(a - b)
        scale = max(abs(a), abs(b), 1.0e-300)
        file_abs = max(file_abs, delta)
        file_rel = max(file_rel, delta / scale)
    overall_abs = max(overall_abs, file_abs)
    overall_rel = max(overall_rel, file_rel)
    print(f"ORACLE file={rel} max_absolute={file_abs:.17g} max_relative={file_rel:.17g}")
    if file_abs > absolute_limit and file_rel > relative_limit:
        raise SystemExit(f"ORACLE_FAIL file={rel}")

print(f"ORACLE_PASS max_absolute={overall_abs:.17g} max_relative={overall_rel:.17g}")
PY
