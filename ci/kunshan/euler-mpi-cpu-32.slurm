#!/bin/bash -l
#SBATCH --job-name=oneflow_euler_mpi_cpu32
#SBATCH --partition=kshcnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=110G
#SBATCH --time=00:20:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -uo pipefail

project="${ONEFLOW_KUNSHAN_PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source_dir="${ONEFLOW_EULER_SOURCE_DIR:-}"
if [ -z "$source_dir" ]; then
  if [ -d "$project/ports/kunshan/oneflow_1d_hip" ]; then
    source_dir="$project/ports/kunshan/oneflow_1d_hip"
  else
    source_dir="$project/src/ports/kunshan/oneflow_1d_hip"
  fi
fi
build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-oneflow-euler}}"
cmake_bin="${CMAKE_BIN:-cmake}"

module purge
module load compiler/devtoolset/7.3.1
module load mpi/hpcx/2.4.1-gcc-7.3.1

export CC="$(command -v gcc)"
export CXX="$(command -v g++)"
mkdir -p "$build_dir"

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DONEFLOW_1D_ENABLE_MPI_BENCHMARK=ON
rc_config=$?
rc_build=1
if [ "$rc_config" -eq 0 ]; then
  "$cmake_bin" --build "$build_dir" --target oneflow_1d_euler_mpi_cpu_benchmark \
    --parallel 32
  rc_build=$?
fi

rc_run=1
if [ "$rc_build" -eq 0 ]; then
  for nx in 65536 262144 1048576 4194304; do
    echo "===== CPU MPI nx=$nx ====="
    mpirun -np 32 --bind-to core \
      "$build_dir/oneflow_1d_euler_mpi_cpu_benchmark" "$nx" 100 2 1
    rc_run=$?
    if [ "$rc_run" -ne 0 ]; then break; fi
  done
fi

echo "RESULT config=$rc_config build=$rc_build run=$rc_run"
if [ "$rc_config" -eq 0 ] && [ "$rc_build" -eq 0 ] && [ "$rc_run" -eq 0 ]; then
  exit 0
fi
exit 1
