#!/bin/bash -l
# Accuracy-gated 3D main-solver CPU/HIP timing on one Kunshan DCU node.
# Reuses an already validated root HIP build; it does not rebuild solver code.
#SBATCH --job-name=oneflow_f3_benchmark
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -uo pipefail

root="$HOME/.scnet-hpc/oneflow"
source_dir="$root/src/OneFLOW-fd8d9eca"
build_dir="$root/builds/root-dcu-fd8d9eca"
artifact_dir="${ONEFLOW_F3_ARTIFACT_DIR:-$root/runs/20260916/f3-main-solver-benchmark}"
runner="$artifact_dir/stage-runner.slurm"
trace_steps="${ONEFLOW_F3_TRACE_STEPS:-1}"
benchmark_repeats="${ONEFLOW_F3_BENCH_REPEATS:-3}"
benchmark_warmup="${ONEFLOW_F3_BENCH_WARMUP:-1}"
cpu_ranks="${ONEFLOW_F3_CPU_RANKS:-1}"
hip_ranks="${ONEFLOW_F3_HIP_RANKS:-1}"

module purge
module load compiler/cmake/3.25.0
module load compiler/gcc/9.3.0
module load mpi/openmpi/gcc-9.3.0/4.1.5
module load compiler/dtk/26.04

cmake_bin="$(command -v cmake)"
ctest_bin="$(dirname "$cmake_bin")/ctest"
clang_bin="$(command -v clang)"
gcc_root="$(dirname "$(dirname "$(command -v gcc)")")"
llvm_root="$(dirname "$(dirname "$clang_bin")")"
dtk_root="$(dirname "$llvm_root")"
amd_comgr_dir="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr"
amd_device_libs_dir="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs"
mpi_root="/public/software/mpi/openmpi/openmpi-4.1.5"
metis_root="$root/deps/metis-install"

export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export CMAKE_PREFIX_PATH="$amd_device_libs_dir:$dtk_root/dcc/comgr/lib64/cmake:$dtk_root:$dtk_root/.hyhal:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:$dtk_root/dcc/lib64:$dtk_root/dcc/comgr/lib64:$metis_root/lib:$mpi_root/lib:${LD_LIBRARY_PATH:-}"
export CC="$llvm_root/bin/clang"
export CXX="$llvm_root/bin/clang++"
export MPI_HOME_INC="$mpi_root/include"
export MPI_HOME_LIB="$mpi_root/lib/libmpi.so"
export METIS_HOME_INC="$metis_root/include"
export METIS_HOME_LIB="$metis_root/lib/libmetis.so"
export CGNS_HOME_INC="/public/software/mathlib/CGNS-4.2.0/src"
export CGNS_HOME_LIB="/public/software/mathlib/CGNS-4.2.0/src/lib/libcgns.a"
export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"

job_tmp="$root/tmp/${SLURM_JOB_ID:-f3-main-solver-benchmark}"
mkdir -p "$job_tmp" "$artifact_dir"
export TMPDIR="$job_tmp"

{
    echo "revision=fd8d9eca"
    echo "date=$(date -Iseconds)"
    echo "partition=${SLURM_JOB_PARTITION:-unset}"
    echo "gres=${SLURM_JOB_GRES:-unset}"
    echo "clang=$("$CXX" --version | head -1)"
    echo "gcc=$(gcc --version | head -1)"
    echo "cmake=$("$cmake_bin" --version | head -1)"
    echo "hip_arch=gfx906"
    echo "hip_visible_devices=${HIP_VISIBLE_DEVICES:-unset}"
    echo "rocr_visible_devices=${ROCR_VISIBLE_DEVICES:-unset}"
    echo "trace_steps=$trace_steps"
    echo "benchmark_repeats=$benchmark_repeats"
    echo "benchmark_warmup=$benchmark_warmup"
    echo "benchmark_cpu_ranks=$cpu_ranks"
    echo "benchmark_hip_ranks=$hip_ranks"
} > "$artifact_dir/environment.txt" 2>&1

if command -v rocminfo >/dev/null 2>&1; then
    rocminfo > "$artifact_dir/rocminfo.txt" 2>&1 || true
fi
if command -v hipconfig >/dev/null 2>&1; then
    hipconfig --full > "$artifact_dir/hipconfig.txt" 2>&1 || true
fi

rc_smoke=1
rc_gtest=1
rc_ctest=1
rc_stage=1

if [ -x "$build_dir/bin/OneFLOWHipSmoke" ]; then
    "$build_dir/bin/OneFLOWHipSmoke" > "$artifact_dir/smoke.log" 2>&1
    rc_smoke=$?
fi

test_binary="$build_dir/bin/oneflow_euler_backend_hip_test"
if [ -x "$test_binary" ]; then
    find "$build_dir" -maxdepth 3 -type f -name "*_tests.cmake" -delete
    "$test_binary" > "$artifact_dir/gtest.log" 2>&1
    rc_gtest=$?
    if [ "$rc_gtest" -eq 0 ] && ! grep -Eq "\[  PASSED  \] [1-9][0-9]* tests?\.?$" "$artifact_dir/gtest.log"; then
        rc_gtest=1
    fi
fi

if [ -x "$ctest_bin" ]; then
    "$source_dir/ci/kunshan/ctest-nonempty.sh" "$build_dir" \
        -L hardware -R HIP \
        --output-on-failure > "$artifact_dir/ctest.log" 2>&1
    rc_ctest=$?
fi

if [ "$rc_smoke" -eq 0 ] && [ "$rc_gtest" -eq 0 ] && [ "$rc_ctest" -eq 0 ] && [ -r "$runner" ]; then
    export ONEFLOW_F2_STAGE_SOURCE_DIR="$source_dir"
    export ONEFLOW_F2_STAGE_BINARY="$build_dir/bin/OneFLOW"
    export ONEFLOW_F2_STAGE_RUN_DIR="$artifact_dir/stage"
    export ONEFLOW_F2_STAGE_CASE="test/cfd/m6wingroe_sa"
    export ONEFLOW_F2_STAGE_STEPS="$trace_steps"
    export ONEFLOW_F2_STAGE_TIME_INTEGRAL=1
    export ONEFLOW_F2_STAGE_RK_STAGES=3
    export ONEFLOW_F2_STAGE_DISABLE_LIMITER=1
    export ONEFLOW_F2_STAGE_CFL=0.01
    export ONEFLOW_F2_STAGE_BENCHMARK=1
    export ONEFLOW_F2_STAGE_BENCH_STEPS="$trace_steps"
    export ONEFLOW_F2_STAGE_BENCH_REPEATS="$benchmark_repeats"
    export ONEFLOW_F2_STAGE_BENCH_WARMUP="$benchmark_warmup"
    export ONEFLOW_F2_STAGE_BENCH_CPU_RANKS="$cpu_ranks"
    export ONEFLOW_F2_STAGE_BENCH_HIP_RANKS="$hip_ranks"
    bash "$runner" > "$artifact_dir/stage-runner.log" 2>&1
    rc_stage=$?
fi

{
    echo "smoke=$rc_smoke"
    echo "hip_gtest=$rc_gtest"
    echo "hip_ctest=$rc_ctest"
    echo "stage_and_benchmark=$rc_stage"
} | tee "$artifact_dir/result.txt"

if [ "$rc_smoke" -eq 0 ] && [ "$rc_gtest" -eq 0 ] && \
   [ "$rc_ctest" -eq 0 ] && [ "$rc_stage" -eq 0 ]; then
    exit 0
fi
exit 1
