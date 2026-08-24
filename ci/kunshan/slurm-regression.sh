#!/bin/bash -l

set -uo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <source-dir> <run-dir> <config-file> <suite>" >&2
    exit 64
fi

source_dir=$1
run_dir=$2
config_file=$3
suite=$4

# shellcheck source=/dev/null
source "$config_file"

if ! declare -F kunshan_load_environment >/dev/null; then
    echo "Remote config must define kunshan_load_environment()." >&2
    exit 78
fi

kunshan_load_environment

build_dir="$run_dir/build"
work_dir="$run_dir/work"
artifact_dir="$run_dir/artifacts"
binary="$build_dir/bin/OneFLOW"

mkdir -p "$build_dir" "$work_dir" "$artifact_dir"

build_jobs=${KUNSHAN_BUILD_JOBS:-16}
serial_timeout=${KUNSHAN_SERIAL_TIMEOUT_SECONDS:-180}
mpi_timeout=${KUNSHAN_MPI_TIMEOUT_SECONDS:-180}

if ! declare -p KUNSHAN_CMAKE_ARGS >/dev/null 2>&1; then
    KUNSHAN_CMAKE_ARGS=()
fi

{
    echo "suite=$suite"
    echo "host=$(hostname)"
    echo "slurm_job_id=${SLURM_JOB_ID:-unset}"
    echo "source_dir=$source_dir"
    echo "run_dir=$run_dir"
    echo "build_jobs=$build_jobs"
} > "$artifact_dir/environment.txt"

cmake -S "$source_dir" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DONEFLOW_ACCEL_BACKEND=CPU \
    -DONEFLOW_ENABLE_MULTI_DEVICE=ON \
    "${KUNSHAN_CMAKE_ARGS[@]}"
config_rc=$?
if [ "$config_rc" -ne 0 ]; then
    echo "CMake configure failed: $config_rc" >&2
    exit "$config_rc"
fi

cmake --build "$build_dir" --parallel "$build_jobs"
build_rc=$?
if [ "$build_rc" -ne 0 ]; then
    echo "Build failed: $build_rc" >&2
    exit "$build_rc"
fi

test -x "$binary" || {
    echo "Missing OneFLOW binary: $binary" >&2
    exit 127
}

rm -rf "$build_dir/bin/system"
cp -a "$source_dir/system" "$build_dir/bin/system"

run_cpu_serial()
{
    local cpu_work="$work_dir/cpu-serial"
    local case_name

    mkdir -p "$cpu_work"
    while IFS= read -r case_name; do
        [ -n "$case_name" ] || continue
        cp -a "$source_dir/test/$case_name" "$cpu_work/$case_name"
        rm -rf \
            "$cpu_work/$case_name/results" \
            "$cpu_work/$case_name/restart" \
            "$cpu_work/$case_name/log"
        mkdir -p \
            "$cpu_work/$case_name/results" \
            "$cpu_work/$case_name/restart" \
            "$cpu_work/$case_name/log"
    done < "$source_dir/test/suites/cpu-serial.txt"

    (
        cd "$cpu_work"
        timeout --signal=TERM --kill-after=10s "$serial_timeout" \
            python3 "$source_dir/test/test.py" \
            "mpirun -np 1" "$binary" "$source_dir/test/suites/cpu-serial.txt"
    )
}

run_mpi4()
{
    : "${KUNSHAN_MPI4_CASE_ROOT:?MPI suite requires KUNSHAN_MPI4_CASE_ROOT}"
    test -d "$KUNSHAN_MPI4_CASE_ROOT"

    local mpi_work="$work_dir/mpi4"
    local case_name=${KUNSHAN_MPI4_CASE_NAME:-m6wingroe_sa_np4}
    local suite_file="$mpi_work/mpi4.txt"

    mkdir -p "$mpi_work"
    cp -a "$KUNSHAN_MPI4_CASE_ROOT" "$mpi_work/$case_name"
    rm -rf \
        "$mpi_work/$case_name/results" \
        "$mpi_work/$case_name/restart" \
        "$mpi_work/$case_name/log"
    mkdir -p \
        "$mpi_work/$case_name/results" \
        "$mpi_work/$case_name/restart" \
        "$mpi_work/$case_name/log"
    printf '%s\n' "$case_name" > "$suite_file"

    (
        cd "$mpi_work"
        export OMP_NUM_THREADS=${KUNSHAN_MPI_OMP_THREADS:-1}
        timeout --signal=TERM --kill-after=10s "$mpi_timeout" \
            python3 "$source_dir/test/test.py" \
            "mpirun -np ${KUNSHAN_MPI_NTASKS:-4}" \
            "$binary" "$suite_file"
    )
}

overall_rc=0
case "$suite" in
    all)
        run_cpu_serial || overall_rc=$?
        if [ "$overall_rc" -eq 0 ]; then
            run_mpi4 || overall_rc=$?
        fi
        ;;
    cpu-serial)
        run_cpu_serial || overall_rc=$?
        ;;
    mpi4)
        run_mpi4 || overall_rc=$?
        ;;
    *)
        echo "Unsupported suite: $suite" >&2
        overall_rc=64
        ;;
esac

{
    echo "suite=$suite"
    echo "result=$overall_rc"
    echo "completed_at=$(date -Is)"
} > "$artifact_dir/summary.txt"

find "$work_dir" -path '*/log/*.log' -type f \
    -exec cp --parents '{}' "$artifact_dir" ';' 2>/dev/null || true

exit "$overall_rc"
