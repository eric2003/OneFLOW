#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <source-dir> <run-dir> <config-file> <suite>" >&2
    exit 64
fi

source_dir=$1
run_dir=$2
config_file=$3
suite=$4

case "$suite" in
    all|cpu-serial|mpi4) ;;
    *)
        echo "Unsupported regression suite: $suite" >&2
        exit 64
        ;;
esac

test -d "$source_dir"
test -d "$run_dir"
test -r "$config_file"

# The remote configuration owns all cluster-specific module and dependency
# paths. It is deliberately not committed to the repository.
# shellcheck source=/dev/null
source "$config_file"

: "${KUNSHAN_CPU_PARTITION:?Missing KUNSHAN_CPU_PARTITION in remote config}"
: "${KUNSHAN_REGRESSION_MEM:?Missing KUNSHAN_REGRESSION_MEM in remote config}"
: "${KUNSHAN_REGRESSION_TIME:?Missing KUNSHAN_REGRESSION_TIME in remote config}"

if [ "$suite" = "cpu-serial" ]; then
    ntasks=${KUNSHAN_SERIAL_NTASKS:-1}
    cpus_per_task=${KUNSHAN_SERIAL_CPUS_PER_TASK:-16}
else
    ntasks=${KUNSHAN_MPI_NTASKS:-4}
    cpus_per_task=${KUNSHAN_MPI_CPUS_PER_TASK:-4}
fi

mkdir -p "$run_dir/artifacts"

job_id=$(sbatch --parsable \
    --job-name=oneflow-regression \
    --partition="$KUNSHAN_CPU_PARTITION" \
    --ntasks="$ntasks" \
    --cpus-per-task="$cpus_per_task" \
    --mem="$KUNSHAN_REGRESSION_MEM" \
    --time="$KUNSHAN_REGRESSION_TIME" \
    --output="$run_dir/artifacts/slurm-%j.out" \
    --error="$run_dir/artifacts/slurm-%j.err" \
    --export=NONE \
    "$source_dir/ci/kunshan/slurm-regression.sh" \
    "$source_dir" "$run_dir" "$config_file" "$suite")

printf '%s\n' "$job_id" > "$run_dir/artifacts/job-id.txt"
echo "Submitted Slurm job: $job_id"

while squeue -h -j "$job_id" | grep -q .; do
    sleep 10
done

sacct -j "$job_id" \
    --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS \
    -P > "$run_dir/artifacts/sacct.txt"

state=$(sacct -n -X -j "$job_id" --format=State | awk 'NF { print $1; exit }')
exit_code=$(sacct -n -X -j "$job_id" --format=ExitCode | awk 'NF { print $1; exit }')

echo "Slurm state: ${state:-UNKNOWN}"
echo "Slurm exit code: ${exit_code:-UNKNOWN}"

if [ "$state" != "COMPLETED" ] || [ "$exit_code" != "0:0" ]; then
    echo "Kunshan regression failed; see downloaded Slurm logs." >&2
    exit 1
fi
