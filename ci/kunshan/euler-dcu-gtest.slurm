#!/bin/bash -l
#SBATCH --job-name=oneflow_euler_dcu_gtest
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
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

build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-oneflow-euler-gtest}}"
artifact_dir="${ONEFLOW_ARTIFACT_DIR:-$project/.oneflow-artifacts/${SLURM_JOB_ID:-local}}"
# Resolve host tools before `module purge`; the DTK module may put an old
# system CMake first on PATH.
cmake_bin="${CMAKE_BIN:-$(command -v cmake)}"
ctest_bin="${CTEST_BIN:-$(dirname "$cmake_bin")/ctest}"
dtk_module="${DTK_MODULE:-compiler/dtk/26.04}"

module purge
module load "$dtk_module"
clang_bin="$(command -v clang)"
llvm_root="${LLVM_ROOT:-$(dirname "$(dirname "$clang_bin")")}"
dtk_root="${DTK_ROOT:-$(dirname "$llvm_root")}"
export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:${LD_LIBRARY_PATH:-}"
export CC="$llvm_root/bin/clang"
export CXX="$llvm_root/bin/clang++"
export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"

job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/${SLURM_JOB_ID:-local}"
mkdir -p "$job_tmp" "$build_dir" "$artifact_dir"
export TMPDIR="$job_tmp"

{
  echo "host=$(hostname)"
  echo "job_id=${SLURM_JOB_ID:-local}"
  echo "source_dir=$source_dir"
  echo "dtk_module=$dtk_module"
  echo "clang=$("$CC" --version | head -1)"
  echo "cmake=$("$cmake_bin" --version | head -1)"
  echo "hip_visible_devices=${HIP_VISIBLE_DEVICES:-unset}"
  echo "rocr_visible_devices=${ROCR_VISIBLE_DEVICES:-unset}"
  echo "slurm_partition=${SLURM_JOB_PARTITION:-unset}"
  echo "slurm_gres=${SLURM_JOB_GRES:-unset}"
} > "$artifact_dir/environment.txt"

if command -v rocminfo >/dev/null 2>&1; then
  rocminfo > "$artifact_dir/rocminfo.txt" 2>&1 || true
fi
if command -v hipconfig >/dev/null 2>&1; then
  hipconfig --full > "$artifact_dir/hipconfig.txt" 2>&1 || true
fi

rc_config=1
rc_build=1
rc_test=1

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_HIP_COMPILER="$CXX" \
  -DCMAKE_PREFIX_PATH="$dtk_root;$dtk_root/.hyhal" \
  -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
  -DHSA_HEADER="$dtk_root/.hyhal/include" \
  -DCMAKE_HIP_ARCHITECTURES="${ONEFLOW_HIP_ARCHITECTURES:-gfx906}" \
  -DONEFLOW_1D_ENABLE_HIP=ON \
  -DONEFLOW_1D_ENABLE_GTEST=ON \
  > "$artifact_dir/configure.log" 2>&1
rc_config=$?

if [ "$rc_config" -eq 0 ]; then
  "$cmake_bin" --build "$build_dir" \
    --target oneflow_1d_euler_cpu_backend_test \
    oneflow_1d_euler_hip_backend_test \
    --parallel 8 > "$artifact_dir/build.log" 2>&1
  rc_build=$?
fi

if [ "$rc_build" -eq 0 ]; then
  # PRE_TEST discovery can retain stale *_tests.cmake files when CMake metadata
  # changes but the test binaries do not. Force regeneration for this run.
  find "$build_dir" -maxdepth 1 -type f -name "*_tests.cmake" -delete
  test_binary="$build_dir/oneflow_1d_euler_hip_backend_test"
  "$test_binary" > "$artifact_dir/gtest.log" 2>&1
  rc_test=$?
  if [ "$rc_test" -eq 0 ] && ! grep -Eq "\[  PASSED  \] 6 tests" "$artifact_dir/gtest.log"; then
    echo "HIP GoogleTest did not report 6 passed tests" >&2
    rc_test=1
  fi
  if [ "$rc_test" -eq 0 ]; then
    "$ctest_bin" --test-dir "$build_dir" \
      -R '^HIP\.' \
      --output-on-failure > "$artifact_dir/ctest.log" 2>&1
    rc_ctest=$?
    if [ "$rc_ctest" -ne 0 ] || ! grep -Eq "100% tests passed, 0 tests failed out of 6" "$artifact_dir/ctest.log"; then
      echo "HIP CTest did not report 6 passed tests" >&2
      rc_test=1
    fi
  fi
fi

{
  echo "config=$rc_config"
  echo "build=$rc_build"
  echo "test=$rc_test"
} | tee "$artifact_dir/result.txt"

if [ "$rc_config" -eq 0 ] && [ "$rc_build" -eq 0 ] && [ "$rc_test" -eq 0 ]; then
  exit 0
fi
exit 1
