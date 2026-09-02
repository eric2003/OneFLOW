#!/bin/bash -l
# 4-rank/4-DCU MPI regression and multi-card benchmark.
# Resource ratio follows the kseshell profile: 4 DCU, 4x8 CPU, 3569 MB/CPU.
#SBATCH --job-name=oneflow_euler_dcu_mpi4_regression
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3569
#SBATCH --gres=dcu:4
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
build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-oneflow-euler}-${SLURM_JOB_ID:-manual}}"
echo "ONEFLOW_PROJECT=$project"
echo "ONEFLOW_SOURCE_DIR=$source_dir"
echo "ONEFLOW_BUILD_DIR=$build_dir"
cmake_bin="${CMAKE_BIN:-cmake}"
dtk_module="${DTK_MODULE:-compiler/dtk/26.04}"
job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/${SLURM_JOB_ID:-manual-$$}"

module purge
module load compiler/devtoolset/7.3.1
module load mpi/hpcx/2.4.1-gcc-7.3.1
module load compiler/cmake/3.25.0
module load "$dtk_module"
clang_bin="$(command -v clang)"
llvm_root="${LLVM_ROOT:-$(dirname "$(dirname "$clang_bin")")}"
dtk_root="${DTK_ROOT:-$(dirname "$llvm_root")}"
export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:${LD_LIBRARY_PATH:-}"
export CMAKE_PREFIX_PATH="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr:$dtk_root/dcc/lib64/cmake/AMDDeviceLibs:$dtk_root:$dtk_root/.hyhal"
hip_compiler="$llvm_root/bin/clang++"
export CC="$(command -v mpicc)"
export CXX="$(command -v mpicxx)"
unset HIP_VISIBLE_DEVICES
unset ROCR_VISIBLE_DEVICES
export ONEFLOW_MPI_REQUIRE_ONE_RANK_PER_DEVICE=1

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
  -DCMAKE_HIP_COMPILER="$hip_compiler" \
  -DCMAKE_PREFIX_PATH="$dtk_root;$dtk_root/.hyhal;$dtk_root/dcc/comgr;$dtk_root/dcc/lib64" \
  -Damd_comgr_DIR="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr" \
  -DAMDDeviceLibs_DIR="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs" \
  -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
  -DHSA_HEADER="$dtk_root/.hyhal/include" \
  -DCMAKE_HIP_ARCHITECTURES=gfx906 \
  -DONEFLOW_1D_ENABLE_HIP=ON \
  -DONEFLOW_1D_ENABLE_MPI_BENCHMARK=ON
)
rc_config=$?
if [ "$rc_config" -eq 0 ]; then
  "$cmake_bin" --build "$build_dir" \
    --target oneflow_1d_euler_mpi_hip_benchmark \
    oneflow_1d_euler_mpi_regression -- -j8
  rc_build=$?
fi

if [ "$rc_build" -eq 0 ]; then
  rc_run=0
  regression_out="$job_tmp/mpi4-regression.out"
  mpirun -np 4 --bind-to core \
    "$build_dir/oneflow_1d_euler_mpi_regression" 65536 20 \
    >"$regression_out" 2>&1
  regression_rc=$?
  cat "$regression_out"
  if [ "$regression_rc" -ne 0 ] \
     || ! grep -q '^result=PASS$' "$regression_out" \
     || ! grep -q 'visible_devices=4' "$regression_out"; then
    echo "DCU_MPI4_REGRESSION=FAIL"
    rc_run=1
  else
    echo "DCU_MPI4_REGRESSION=PASS"
  fi

  if [ "$rc_run" -eq 0 ]; then
    for nx in 65536 262144 1048576 4194304; do
      benchmark_out="$job_tmp/mpi4-benchmark-${nx}.out"
      echo "===== DCU MPI4 benchmark nx=$nx ====="
      mpirun -np 4 --bind-to core \
        "$build_dir/oneflow_1d_euler_mpi_hip_benchmark" "$nx" 100 2 1 \
        >"$benchmark_out" 2>&1
      benchmark_rc=$?
      cat "$benchmark_out"
      if [ "$benchmark_rc" -ne 0 ] \
         || ! grep -q 'visible_devices=4' "$benchmark_out"; then
        echo "DCU_MPI4_BENCHMARK=FAIL nx=$nx"
        rc_run=1
        break
      fi
      echo "DCU_MPI4_BENCHMARK=PASS nx=$nx"
    done
  fi
fi

echo "RESULT config=$rc_config build=$rc_build run=$rc_run"
if [ "$rc_config" -eq 0 ] && [ "$rc_build" -eq 0 ] && [ "$rc_run" -eq 0 ]; then
  exit 0
fi
exit 1
