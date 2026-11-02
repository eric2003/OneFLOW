#!/bin/bash -l
# Build and contract-test the root 3D HIP solver on a Kunshan DCU node.
#
# The root C/C++ project uses the current GCC 16 toolchain.  DTK 26.04
# clang 17 currently needs GCC 9 libstdc++ headers for HIP compiler probing;
# the final link is explicitly directed to the GCC 16 runtime so objects and
# third-party libraries use one ABI/runtime family.
#SBATCH --job-name=oneflow_f2_build
#SBATCH --partition=kshdnormal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=27G
#SBATCH --gres=dcu:1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

: "${ONEFLOW_F2_BUILD_SOURCE_DIR:?Set ONEFLOW_F2_BUILD_SOURCE_DIR}"
: "${ONEFLOW_F2_BUILD_DIR:?Set ONEFLOW_F2_BUILD_DIR}"

module purge
module load gcc/16.2.0
module load openmpi/5.0.11
module load cmake/4.4.3
module load compiler/dtk/26.04
module load python/3.8.10

export OMPI_MCA_coll="${OMPI_MCA_coll:-^hcoll}"

clang_bin="$(command -v clang)"
llvm_root="$(dirname "$(dirname "$clang_bin")")"
dtk_root="$(dirname "$llvm_root")"
mpi_root="$(dirname "$(dirname "$(command -v mpirun)")")"
gcc16_root="$(dirname "$(dirname "$(command -v gcc)")")"
gcc16_lib="$gcc16_root/lib64"
gcc9_root="${ONEFLOW_HIP_HOST_GCC_ROOT:-/public/software/compiler/gcc-9.3.0}"
metis_root="${ONEFLOW_F2_BUILD_METIS_ROOT:-$HOME/.scnet-hpc/oneflow/deps/metis-install}"
amd_comgr_dir="$dtk_root/dcc/comgr/lib64/cmake/amd_comgr"
amd_device_libs_dir="$dtk_root/dcc/lib64/cmake/AMDDeviceLibs"

test -d "$gcc9_root"
test -d "$gcc16_lib"
test -d "$metis_root"

export PATH="$dtk_root/bin:$llvm_root/bin:$PATH"
export CMAKE_PREFIX_PATH="$amd_device_libs_dir:$dtk_root/dcc/comgr/lib64/cmake:$dtk_root:$dtk_root/.hyhal${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export LD_LIBRARY_PATH="$gcc16_lib:$dtk_root/hip/lib:$llvm_root/lib:$llvm_root/lib64:$dtk_root/lib:$dtk_root/lib64:$dtk_root/dcc/lib64:$dtk_root/dcc/comgr/lib64:$metis_root/lib:$mpi_root/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export HIP_VISIBLE_DEVICES="${HIP_VISIBLE_DEVICES:-0}"
export ROCR_VISIBLE_DEVICES="${ROCR_VISIBLE_DEVICES:-0}"

export MPI_HOME_INC="$mpi_root/include"
export MPI_HOME_LIB="$mpi_root/lib/libmpi.so"
export METIS_HOME_INC="$metis_root/include"
export METIS_HOME_LIB="$metis_root/lib/libmetis.so"
export CGNS_HOME_INC="${CGNS_HOME_INC:-/public/software/mathlib/CGNS-4.2.0/src}"
export CGNS_HOME_LIB="${CGNS_HOME_LIB:-/public/software/mathlib/CGNS-4.2.0/src/lib/libcgns.a}"

export CC="${CC:-gcc}"
export CXX="${CXX:-g++}"

link_flags="-L$gcc16_lib -Wl,-rpath,$gcc16_lib -Wl,-rpath-link,$gcc16_lib"
build_log="$ONEFLOW_F2_BUILD_DIR/build.log"
mkdir -p "$ONEFLOW_F2_BUILD_DIR"

# Fail before CMake when the allocation has no usable DCU runtime.  This
# keeps /dev/kfd, architecture, and amd_comgr problems distinct from build
# or solver failures and leaves a bounded diagnostic artifact.
export AMD_COMGR_DIR="$amd_comgr_dir"
"$ONEFLOW_F2_BUILD_SOURCE_DIR/ci/kunshan/dcu-preflight.sh" \
    "$ONEFLOW_F2_BUILD_DIR/dcu-preflight.txt"

{
    echo "source=$ONEFLOW_F2_BUILD_SOURCE_DIR"
    echo "build=$ONEFLOW_F2_BUILD_DIR"
    echo "gcc=$(gcc --version | head -1)"
    echo "hip_clang=$("$llvm_root/bin/clang++" --version | head -1)"
    echo "hip_host_headers=$gcc9_root"
    echo "mpi=$(mpirun --version | head -1)"
    echo "cmake=$(cmake --version | head -1)"
    echo "mpirun=$(command -v mpirun)"
    echo "mpi_root=$mpi_root"
    echo "gcc16_lib=$gcc16_lib"
    echo "OMPI_MCA_coll=$OMPI_MCA_coll"
    echo "hip_arch=${ONEFLOW_HIP_ARCHITECTURES:-gfx906}"
    echo "hip_visible_devices=$HIP_VISIBLE_DEVICES"
    echo "rocr_visible_devices=$ROCR_VISIBLE_DEVICES"
} | tee "$ONEFLOW_F2_BUILD_DIR/environment.txt"

cmake -S "$ONEFLOW_F2_BUILD_SOURCE_DIR" -B "$ONEFLOW_F2_BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="$CC" \
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_HIP_COMPILER="$llvm_root/bin/clang++" \
    -DCMAKE_HIP_FLAGS="--gcc-toolchain=$gcc9_root -stdlib=libstdc++ -L$gcc16_lib" \
    -DCMAKE_EXE_LINKER_FLAGS="$link_flags" \
    -DCMAKE_SHARED_LINKER_FLAGS="$link_flags" \
    -DCMAKE_HIP_ARCHITECTURES="${ONEFLOW_HIP_ARCHITECTURES:-gfx906}" \
    -DCMAKE_PREFIX_PATH="$amd_device_libs_dir;$dtk_root/dcc/comgr/lib64/cmake;$dtk_root;$dtk_root/.hyhal" \
    -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
    -DCMAKE_TRY_COMPILE_PLATFORM_VARIABLES="CMAKE_PREFIX_PATH;amd_comgr_DIR;AMDDeviceLibs_DIR" \
    -Damd_comgr_DIR="$amd_comgr_dir" \
    -DAMDDeviceLibs_DIR="$amd_device_libs_dir" \
    -Dhsa-runtime64_DIR="$dtk_root/.hyhal/lib/cmake/hsa-runtime64" \
    -DHSA_HEADER="$dtk_root/.hyhal/include" \
    -DONEFLOW_ENABLE_HIP=ON \
    -DONEFLOW_ENABLE_HIP_TESTS=ON \
    -DMPI_ENABLE=ON \
    -DMETIS_ENABLE=ON \
    -DCGNS_ENABLE=ON

cmake --build "$ONEFLOW_F2_BUILD_DIR" \
    --target OneFLOW OneFLOWHipSmoke oneflow_euler_backend_hip_test \
             oneflow_hip_resident_state_contract_test \
    --parallel "${ONEFLOW_F2_BUILD_PARALLEL:-8}" 2>&1 | tee "$build_log"

# The solver resolves its runtime "system" directory relative to the
# executable. A targeted build does not install it, so place it explicitly
# instead of leaving a build that cannot start a case.
install -d "$ONEFLOW_F2_BUILD_DIR/bin"
if [ ! -r "$ONEFLOW_F2_BUILD_DIR/bin/system/action/actionFileList.txt" ]; then
    cp -a "$ONEFLOW_F2_BUILD_SOURCE_DIR/system" \
        "$ONEFLOW_F2_BUILD_DIR/bin/system"
fi
test -r "$ONEFLOW_F2_BUILD_DIR/bin/system/action/actionFileList.txt"

"$ONEFLOW_F2_BUILD_DIR/bin/OneFLOWHipSmoke"
"$ONEFLOW_F2_BUILD_DIR/bin/oneflow_euler_backend_hip_test"
# Backend-neutral resident-state generation contract (runs on the CPU).
"$ONEFLOW_F2_BUILD_DIR/bin/oneflow_hip_resident_state_contract_test"

"$ONEFLOW_F2_BUILD_SOURCE_DIR/ci/kunshan/ctest-nonempty.sh" \
    "$ONEFLOW_F2_BUILD_DIR" -L hardware -R HIP \
    --output-on-failure | tee "$ONEFLOW_F2_BUILD_DIR/ctest.log"
