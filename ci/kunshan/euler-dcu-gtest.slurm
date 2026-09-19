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
dtk_module="${DTK_MODULE:-compiler/dtk/26.04}"
gcc_module="${GCC_MODULE:-compiler/gcc/9.3.0}"
cmake_module="${CMAKE_MODULE:-compiler/cmake/3.25.0}"

module purge
module load "$cmake_module"
module load "$gcc_module"
module load "$dtk_module"
# Resolve CMake after module setup. The login environment exposes CMake 2.8,
# which does not support the -S/-B configure form used below.
cmake_bin="${CMAKE_BIN:-$(command -v cmake)}"
ctest_bin="${CTEST_BIN:-$(dirname "$cmake_bin")/ctest}"
clang_bin="$(command -v clang)"
gcc_root="${GCC_ROOT:-$(dirname "$(dirname "$(command -v gcc)")")}"
gcc_toolchain_flag="--gcc-toolchain=$gcc_root"
stdlib_flag="-stdlib=libstdc++"
compiler_flags="$gcc_toolchain_flag $stdlib_flag"
llvm_root="${LLVM_ROOT:-$(dirname "$(dirname "$clang_bin")")}"
dtk_root="${DTK_ROOT:-$(dirname "$llvm_root")}"
amd_comgr_dir="${AMD_COMGR_DIR:-$dtk_root/dcc/comgr/lib64/cmake/amd_comgr}"
amd_device_libs_dir="${AMD_DEVICE_LIBS_DIR:-$dtk_root/dcc/lib64/cmake/AMDDeviceLibs}"
export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export CMAKE_PREFIX_PATH="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs:$dtk_root/dcc/comgr/lib64/cmake:$dtk_root:$dtk_root/.hyhal:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:$dtk_root/dcc/lib64:$dtk_root/dcc/comgr/lib64:${LD_LIBRARY_PATH:-}"
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
  -DCMAKE_CXX_FLAGS="$compiler_flags" \
  -DCMAKE_HIP_FLAGS="$compiler_flags" \
  -DCMAKE_HIP_COMPILER="$CXX" \
  -DCMAKE_PREFIX_PATH="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs;$dtk_root/dcc/comgr/lib64/cmake;$dtk_root;$dtk_root/.hyhal" \
  -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
  -DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES="CMAKE_PREFIX_PATH;amd_comgr_DIR;AMDDeviceLibs_DIR" \
  -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
  -DHSA_HEADER="$dtk_root/.hyhal/include" \
  -Damd_comgr_DIR="$amd_comgr_dir" \
  -DAMDDeviceLibs_DIR="$amd_device_libs_dir" \
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
  if [ "$rc_test" -eq 0 ] && ! grep -Eq "\[  PASSED  \] [1-9][0-9]* tests?\.?$" "$artifact_dir/gtest.log"; then
    echo "HIP GoogleTest did not report a passing test suite" >&2
    rc_test=1
  fi
  if [ "$rc_test" -eq 0 ]; then
    "$ctest_bin" --test-dir "$build_dir" \
      -L 'hardware' \
      -R 'HIP' \
      --output-on-failure > "$artifact_dir/ctest.log" 2>&1
    rc_ctest=$?
    if [ "$rc_ctest" -ne 0 ] || ! grep -Eq "100% tests passed, 0 tests failed out of [1-9][0-9]*" "$artifact_dir/ctest.log"; then
      echo "HIP CTest did not report a fully passing test suite" >&2
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
