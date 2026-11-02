#!/bin/bash -l
# Single-node CPU-rank/DCU-rank scaling matrix for the standalone Euler MPI port.
# Every configuration checks active MPI ranks, decomposition, halo counts and hash.
# HIP runs use the DTK compatibility toolchain; this is not the CPU/MPI default.
#SBATCH --job-name=oneflow_euler_mpi_matrix
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3569
#SBATCH --gres=dcu:4
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

project="${ONEFLOW_KUNSHAN_PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source_dir="${ONEFLOW_EULER_SOURCE_DIR:-$project/ports/kunshan/oneflow_1d_hip}"
build_dir="${ONEFLOW_BUILD_DIR:-$project/.oneflow-build/${SLURM_JOB_NAME:-matrix}-${SLURM_JOB_ID:-manual}}"
job_tmp="${ONEFLOW_JOB_TMP_ROOT:-${HOME:-$project/.tmp}}/.scnet-hpc/tmp/${SLURM_JOB_ID:-manual-$$}"

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
    oneflow_1d_euler_mpi_hip_benchmark -- -j8

run_cpu() {
    local ranks="$1" nx="$2" output="$3"
    mpirun -np "$ranks" --map-by "ppr:${ranks}:node" --bind-to core \
        "$build_dir/oneflow_1d_euler_mpi_cpu_benchmark" "$nx" 100 2 1 \
        > "$output" 2>&1
}

run_hip() {
    local ranks="$1" nx="$2" output="$3"
    case "$ranks" in
        1)
            export HIP_VISIBLE_DEVICES=0
            export ROCR_VISIBLE_DEVICES=0
            ;;
        2)
            export HIP_VISIBLE_DEVICES=0,1
            export ROCR_VISIBLE_DEVICES=0,1
            ;;
        4)
            unset HIP_VISIBLE_DEVICES ROCR_VISIBLE_DEVICES
            ;;
        *)
            echo "unsupported HIP rank/DCU count: $ranks" >&2
            return 2
            ;;
    esac
    mpirun -np "$ranks" --map-by "ppr:${ranks}:node" --bind-to core \
        "$build_dir/oneflow_1d_euler_mpi_hip_benchmark" "$nx" 100 2 1 \
        > "$output" 2>&1
}

check_output() {
    local output="$1" backend="$2" ranks="$3" nx="$4"
    grep -q "^backend=${backend} ranks=${ranks} local_ranks=${ranks} " "$output"
    grep -q '^decomposition=PASS$' "$output"
    grep -q "decomposition_active_ranks=${ranks} " "$output"
    grep -q "global_cells_sum=${nx} " "$output"
    grep -Eq 'local_nx_min=[1-9][0-9]* local_nx_max=[1-9][0-9]*' "$output"
    grep -Eq 'halo_exchanges_sum=[1-9][0-9]* expected_halo_exchanges=[1-9][0-9]*' \
        "$output"
    if [ "$backend" = "HIP_MPI_SHARED_DEVICE" ]; then
        grep -q \
            "^device=.* local_rank=0 local_ranks=${ranks} device_index=0 visible_devices=${ranks}$" \
            "$output"
    fi
}

for nx in 65536 262144 1048576 4194304; do
    for cpu_ranks in 8 16 32; do
        cpu_output="$job_tmp/cpu-${cpu_ranks}-${nx}.out"
        run_cpu "$cpu_ranks" "$nx" "$cpu_output"
        check_output "$cpu_output" CPU_MPI "$cpu_ranks" "$nx"
        cpu_hash="$(sed -n 's/^final_hash=//p' "$cpu_output")"
        cpu_ms="$(sed -n 's/^lifecycle_max_ms=\([0-9.]*\).*/\1/p' "$cpu_output")"
        for hip_ranks in 1 2 4; do
            hip_output="$job_tmp/hip-${hip_ranks}-${nx}.out"
            run_hip "$hip_ranks" "$nx" "$hip_output"
            check_output "$hip_output" HIP_MPI_SHARED_DEVICE "$hip_ranks" "$nx"
            hip_hash="$(sed -n 's/^final_hash=//p' "$hip_output")"
            test "$cpu_hash" = "$hip_hash"
            hip_ms="$(sed -n 's/^lifecycle_max_ms=\([0-9.]*\).*/\1/p' "$hip_output")"
            ratio="$(awk -v c="$cpu_ms" -v h="$hip_ms" \
                'BEGIN { printf "%.6f", c / h }')"
            echo "MATRIX_PASS nodes=1 nx=$nx cpu_ranks=$cpu_ranks \
cpu_cores=$cpu_ranks hip_ranks=$hip_ranks hip_dcus=$hip_ranks \
cpu_ms=$cpu_ms hip_ms=$hip_ms cpu_over_hip=${ratio}x \
active_cpu_ranks=$cpu_ranks active_hip_ranks=$hip_ranks \
cpu_local_ranks=$cpu_ranks hip_local_ranks=$hip_ranks \
hash_match=1 decomposition=PASS"
        done
    done
done
