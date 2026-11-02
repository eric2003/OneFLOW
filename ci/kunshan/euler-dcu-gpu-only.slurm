#!/bin/bash -l
#SBATCH --job-name=oneflow_euler_dcu_gpu_only
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
#SBATCH --time=00:15:00
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
export HIP_VISIBLE_DEVICES=0
export ROCR_VISIBLE_DEVICES=0

job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/$SLURM_JOB_ID"
mkdir -p "$job_tmp" "$build_dir"
export TMPDIR="$job_tmp"

rc_config=1
rc_build=1
rc_correctness=1
rc_gpu_only=1

"$cmake_bin" -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_HIP_COMPILER="$CXX" \
  -DCMAKE_PREFIX_PATH="$dtk_root;$dtk_root/.hyhal" \
  -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
  -DHSA_HEADER="$dtk_root/.hyhal/include" \
  -DCMAKE_HIP_ARCHITECTURES=gfx906 \
  -DONEFLOW_1D_ENABLE_HIP=ON
rc_config=$?

if [ "$rc_config" -eq 0 ]; then
  "$cmake_bin" --build "$build_dir" \
    --target oneflow_1d_euler_hip oneflow_1d_euler_benchmark \
    oneflow_1d_euler_stateful_benchmark oneflow_1d_euler_gpu_only_benchmark --parallel 8
  rc_build=$?
fi

if [ "$rc_build" -eq 0 ]; then
  "$build_dir/oneflow_1d_euler_hip"
  rc_correctness=$?
  if [ "$rc_correctness" -eq 0 ]; then
    echo "===== OneFLOW GPU-only sustained benchmark ====="
    rocm-smi -d 0 -sdmon u -ddmon 1 -cdmon 15 > "$job_tmp/rocm-smi.csv" 2>&1 &
    smi_pid=$!
    "$build_dir/oneflow_1d_euler_gpu_only_benchmark" 4194304 100 20 2
    rc_gpu_only=$?
    kill "$smi_pid" 2>/dev/null || true
    wait "$smi_pid" 2>/dev/null || true
  fi
fi

if [ "$rc_config" -eq 0 ] && [ "$rc_build" -eq 0 ] && [ "$rc_correctness" -eq 0 ] && [ "$rc_gpu_only" -eq 0 ]; then
  echo "RESULT config=$rc_config build=$rc_build correctness=$rc_correctness gpu_only=$rc_gpu_only"
  exit 0
fi

echo "RESULT config=$rc_config build=$rc_build correctness=$rc_correctness gpu_only=$rc_gpu_only"
exit 1
