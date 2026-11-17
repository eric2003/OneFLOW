#!/bin/bash -l
# Two-node device-pointer MPI probe for Hygon DCU/HIP.
# This is a capability gate, not a solver benchmark.
# It must pass before ONEFLOW_GPU_AWARE_MPI_CONFIRMED=1 is used.
#SBATCH --job-name=oneflow_gpu_aware_mpi_probe
#SBATCH --partition=kshdnormal
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=dcu:2
#SBATCH --time=00:15:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

project="${ONEFLOW_KUNSHAN_PROJECT:-${SLURM_SUBMIT_DIR:-$PWD}}"
source_file="${ONEFLOW_GPU_AWARE_PROBE_SOURCE:-$project/ci/kunshan/gpu-aware-mpi-probe.cpp}"
build_dir="${ONEFLOW_GPU_AWARE_PROBE_BUILD:-$project/.oneflow-build/gpu-aware-mpi-probe-${SLURM_JOB_ID:-manual}}"

if [[ ! -f "$source_file" ]]; then
  echo "GPU-aware MPI probe source not found: $source_file" >&2
  exit 2
fi

module purge
module load gcc/16.2.0
module load openmpi/5.0.11
module load cmake/4.4.3
module load compiler/dtk/26.04

export OMPI_MCA_coll=^hcoll
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
unset OMPI_MCA_opal_cuda_support
unset OMPI_MCA_mpi_cuda_support

mkdir -p "$build_dir"
probe_bin="$build_dir/gpu-aware-mpi-probe"

# hipcc supplies the HIP platform/compiler definitions; mpicxx contributes
# only the selected MPI include and link flags.
hipcc -std=c++17 -O2 "$source_file" \
  -o "$probe_bin" \
  $(mpicxx --showme:compile) \
  $(mpicxx --showme:link)

mpirun -np 2 \
  --map-by ppr:1:node \
  --bind-to core \
  "$probe_bin"

echo "GPU_AWARE_MPI_TWO_NODE_PROBE_PASS"
