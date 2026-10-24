#!/bin/bash -l
# Two-node, four-rank HIP MPI regression and benchmark.
# This is a HIP-specific DTK compatibility runner. It is not the CPU/MPI
# default toolchain used by the standard CPU suites.
#SBATCH --job-name=oneflow_euler_dcu_mpi4_2node
#SBATCH --partition=kshdnormal
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3569
#SBATCH --gres=dcu:2
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -o pipefail

project="${ONEFLOW_KUNSHAN_PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source_dir="${ONEFLOW_EULER_SOURCE_DIR:-$project/ports/kunshan/oneflow_1d_hip}"
build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-oneflow-euler}-${SLURM_JOB_ID:-manual}}"
job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/${SLURM_JOB_ID:-manual-$$}"
mpi_args=( -np 4 --map-by ppr:2:node --bind-to core )

if [ ! -f "$source_dir/CMakeLists.txt" ]; then
    echo "Euler source CMakeLists.txt not found: $source_dir" >&2
    exit 2
fi

module purge
module load compiler/devtoolset/7.3.1
module load mpi/hpcx/2.4.1-gcc-7.3.1
module load compiler/cmake/3.25.0
module load compiler/dtk/26.04

clang_bin="$(command -v clang)"
llvm_root="$(dirname "$(dirname "$clang_bin")")"
dtk_root="$(dirname "$llvm_root")"
export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export LD_LIBRARY_PATH="$dtk_root/hip/lib:$llvm_root/lib:$dtk_root/lib:$dtk_root/lib64:${LD_LIBRARY_PATH:-}"
export CMAKE_PREFIX_PATH="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr:$dtk_root/dcc/lib64/cmake/AMDDeviceLibs:$dtk_root:$dtk_root/.hyhal"
export OMPI_MCA_coll="${OMPI_MCA_coll:-^hcoll}"
export ONEFLOW_MPI_REQUIRE_ONE_RANK_PER_DEVICE=1
export ONEFLOW_MPI_GPU_AWARE=1
export ONEFLOW_MPI_GPU_AWARE_CONFIRMED=1

mkdir -p "$build_dir" "$job_tmp"
export TMPDIR="$job_tmp"

cmake -S "$source_dir" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$(command -v mpicc)" \
    -DCMAKE_CXX_COMPILER="$(command -v mpicxx)" \
    -DCMAKE_HIP_COMPILER="$llvm_root/bin/clang++" \
    -DCMAKE_PREFIX_PATH="$dtk_root;$dtk_root/.hyhal;$dtk_root/dcc/comgr;$dtk_root/dcc/lib64" \
    -Damd_comgr_DIR="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr" \
    -DAMDDeviceLibs_DIR="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs" \
    -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
    -DHSA_HEADER="$dtk_root/.hyhal/include" \
    -DCMAKE_HIP_ARCHITECTURES=gfx906 \
    -DONEFLOW_1D_ENABLE_HIP=ON \
    -DONEFLOW_1D_ENABLE_MPI_BENCHMARK=ON
cmake --build "$build_dir" \
    --target oneflow_1d_euler_mpi_cpu_benchmark \
    oneflow_1d_euler_mpi_hip_benchmark \
    oneflow_1d_euler_mpi_regression -- -j8

regression_out="$job_tmp/regression.out"
mpirun "${mpi_args[@]}" \
    "$build_dir/oneflow_1d_euler_mpi_regression" 65536 20 \
    > "$regression_out" 2>&1
cat "$regression_out"
grep -q '^result=PASS$' "$regression_out"
grep -q '^decomposition=PASS$' "$regression_out"
grep -q 'gpu_aware_mpi=1' "$regression_out"

for nx in 65536 262144 1048576 4194304; do
    cpu_out="$job_tmp/cpu-${nx}.out"
    hip_out="$job_tmp/hip-${nx}.out"
    mpirun "${mpi_args[@]}" \
        "$build_dir/oneflow_1d_euler_mpi_cpu_benchmark" "$nx" 100 2 1 \
        > "$cpu_out" 2>&1
    mpirun "${mpi_args[@]}" \
        "$build_dir/oneflow_1d_euler_mpi_hip_benchmark" "$nx" 100 2 1 \
        > "$hip_out" 2>&1
    cat "$cpu_out"
    cat "$hip_out"
    grep -q '^decomposition=PASS$' "$cpu_out"
    grep -q '^decomposition=PASS$' "$hip_out"
    cpu_hash="$(sed -n 's/^final_hash=//p' "$cpu_out")"
    hip_hash="$(sed -n 's/^final_hash=//p' "$hip_out")"
    test "$cpu_hash" = "$hip_hash"
    echo "DCU_MPI4_2NODE_BENCHMARK=PASS nx=$nx gpu_aware=1"
done

echo "DCU_MPI4_2NODE_REGRESSION=PASS"
