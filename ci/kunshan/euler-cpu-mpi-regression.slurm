#!/bin/bash -l
# 32-rank CPU MPI regression: compare one-rank and 32-rank global hashes.
# Resource ratio follows the kseshell profile: 32 CPU cores, 3569 MB/CPU.
#SBATCH --job-name=oneflow_euler_cpu_mpi_regression
#SBATCH --partition=kshcnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3569
#SBATCH --time=00:40:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -o pipefail

project="${ONEFLOW_KUNSHAN_PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source_dir="${ONEFLOW_EULER_SOURCE_DIR:-}"
if [ -n "$source_dir" ] && [ ! -f "$source_dir/CMakeLists.txt" ]; then
  source_dir=""
fi
if [ -z "$source_dir" ]; then
  if [ -f "$project/ports/kunshan/oneflow_1d_hip/CMakeLists.txt" ]; then
    source_dir="$project/ports/kunshan/oneflow_1d_hip"
  else
    source_dir="$project/src/ports/kunshan/oneflow_1d_hip"
  fi
fi
if [ ! -f "$source_dir/CMakeLists.txt" ]; then
  echo "Euler source CMakeLists.txt not found: $source_dir" >&2
  exit 2
fi
build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-oneflow-euler}}"
echo "ONEFLOW_PROJECT=$project"
echo "ONEFLOW_SOURCE_DIR=$source_dir"
echo "ONEFLOW_BUILD_DIR=$build_dir"
cmake_bin="${CMAKE_BIN:-cmake}"
job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/${SLURM_JOB_ID:-manual-$$}"

module purge
module load compiler/cmake/3.25.0
module load compiler/devtoolset/7.3.1
module load mpi/hpcx/2.4.1-gcc-7.3.1

export CC="$(command -v gcc)"
export CXX="$(command -v g++)"
mkdir -p "$build_dir" "$job_tmp"
export TMPDIR="$job_tmp"

rc_config=1
rc_build=1
rc_run=1

(
  cd "$build_dir" && "$cmake_bin" "$source_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DONEFLOW_1D_ENABLE_MPI_BENCHMARK=ON
)
rc_config=$?
if [ "$rc_config" -eq 0 ]; then
  "$cmake_bin" --build "$build_dir" \
    --target oneflow_1d_euler_mpi_cpu_benchmark -- -j32
  rc_build=$?
fi

if [ "$rc_build" -eq 0 ]; then
  rc_run=0
  for nx in 65536 262144 1048576 4194304; do
    one_rank="$job_tmp/cpu-${nx}-1.out"
    many_rank="$job_tmp/cpu-${nx}-32.out"
    echo "===== CPU MPI regression nx=$nx ====="
    mpirun -np 1 --bind-to core \
      "$build_dir/oneflow_1d_euler_mpi_cpu_benchmark" "$nx" 100 1 1 \
      >"$one_rank" 2>&1
    one_rc=$?
    cat "$one_rank"
    mpirun -np 32 --bind-to core \
      "$build_dir/oneflow_1d_euler_mpi_cpu_benchmark" "$nx" 100 1 1 \
      >"$many_rank" 2>&1
    many_rc=$?
    cat "$many_rank"
    one_hash="$(sed -n 's/^final_hash=\(0x[0-9a-fA-F]*\).*/\1/p' "$one_rank")"
    many_hash="$(sed -n 's/^final_hash=\(0x[0-9a-fA-F]*\).*/\1/p' "$many_rank")"
    echo "hash_1_rank=$one_hash hash_32_rank=$many_hash"
    if [ "$one_rc" -ne 0 ] || [ "$many_rc" -ne 0 ] \
       || [ -z "$one_hash" ] || [ "$one_hash" != "$many_hash" ]; then
      echo "CPU_MPI_REGRESSION=FAIL nx=$nx"
      rc_run=1
      break
    fi
    echo "CPU_MPI_REGRESSION=PASS nx=$nx"
  done
fi

echo "RESULT config=$rc_config build=$rc_build run=$rc_run"
if [ "$rc_config" -eq 0 ] && [ "$rc_build" -eq 0 ] && [ "$rc_run" -eq 0 ]; then
  exit 0
fi
exit 1
