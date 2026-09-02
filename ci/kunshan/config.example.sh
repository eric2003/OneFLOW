#!/usr/bin/env bash

# Copy this file to a private, readable path on the Kunshan login filesystem.
# Do not commit the real configuration: compiler, dependency, fixture, and
# account paths are cluster-specific.

KUNSHAN_CPU_PARTITION="replace-with-cpu-partition"
KUNSHAN_GRES=""
KUNSHAN_REGRESSION_MEM="replace-with-valid-memory-request"
KUNSHAN_REGRESSION_TIME="00:15:00"

KUNSHAN_SERIAL_NTASKS=1
KUNSHAN_SERIAL_CPUS_PER_TASK=16
KUNSHAN_MPI_NTASKS=4
KUNSHAN_MPI_CPUS_PER_TASK=4
KUNSHAN_MPI_OMP_THREADS=1

# CPU is the current scheduled profile. External accelerator profiles may set
# HIP, CUDA, or KOKKOS after providing the matching modules and Slurm resources.
KUNSHAN_ACCEL_BACKEND="CPU"
KUNSHAN_BUILD_JOBS=16
KUNSHAN_SERIAL_TIMEOUT_SECONDS=180
KUNSHAN_MPI_TIMEOUT_SECONDS=180

# This private fixture contains the partitioned grid and its independent
# four-rank CPU baseline. It must be a complete OneFLOW test case directory.
KUNSHAN_MPI4_CASE_ROOT="/dedicated/oneflow-ci/fixtures/m6wingroe_sa_np4"
KUNSHAN_MPI4_CASE_NAME="m6wingroe_sa_np4"

kunshan_load_environment()
{
    module purge
    module load replace-with-compiler-module
    module load replace-with-mpi-module
    module load replace-with-cmake-module
    module load replace-with-python-module

    # Export dependency library paths here when they are not supplied by
    # modules. Do not print credentials or private paths in this function.
}

KUNSHAN_CMAKE_ARGS=(
    -DMPI_ENABLE=ON
    -DMETIS_ENABLE=ON
    -DCGNS_ENABLE=ON
    # Add compiler and dependency include/library cache entries here.
)
