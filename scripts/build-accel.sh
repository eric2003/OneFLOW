#!/usr/bin/env bash
#
# Optional: configure and build OneFLOW with an accelerator SDK (HIP/DCU or CUDA).
#
# This is the build-time half of the accelerator interface:
#
#   runtime backend   ONEFLOW_ACCEL_BACKEND = AUTO | CPU | HIP | CUDA | KOKKOS
#   build-time SDK    ONEFLOW_ACCEL_SDK_ROOT (+ optional _ARCH / _FLAGS), kind inferred
#
# Plain CPU/CI builds never need this script; without an SDK root the CMake
# resolver does nothing. With a root, cmake/OneFLOWAccelSDK.cmake derives
# everything the vendor package expects (for HIP/DCU: ROCM_PATH, the
# amd_comgr / AMDDeviceLibs / hsa-runtime64 directories, and the
# CMAKE_TRY_COMPILE_PLATFORM_VARIABLES that enable_language(HIP) needs), so a
# DCU build is one command instead of a rediscovered set of variables.
#
# The SDK kind is inferred (ONEFLOW_ENABLE_HIP/CUDA, or the prefix contents), so
# there is no kind argument to get wrong.
#
# Usage:
#   scripts/build-accel.sh [--sdk-root <prefix>]
#                          [--arch <gfx906|sm_80>] [--flags "<compiler flags>"]
#                          [--no-sdk-tests] [--jobs N] [--build-dir <dir>]
#                          [--configure-only]
#
# --sdk-root defaults to $ONEFLOW_ACCEL_SDK_ROOT, then to the vendor variable
# ($ROCM_PATH for hip, $CUDAToolkit_ROOT/$CUDA_HOME for cuda). Run it where the
# accelerator is visible when you want CMake's device detection to choose the
# architecture.

set -Eeuo pipefail

# Path-based on purpose: the script must work inside an extracted source
# snapshot (no .git) and with older git builds that lack `git -C`.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "$repo_root"

sdk_root="${ONEFLOW_ACCEL_SDK_ROOT:-}"
arch=""
sdk_flags="${ONEFLOW_ACCEL_SDK_FLAGS:-}"
sdk_tests="ON"
jobs=""
build_dir="build"
configure_only=false

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

usage() {
    sed -n "2,30p" "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --sdk-root) sdk_root="${2:?--sdk-root needs a path}"; shift 2 ;;
        --arch) arch="${2:?--arch needs a value}"; shift 2 ;;
        --flags) sdk_flags="${2:?--flags needs a value}"; shift 2 ;;
        --no-sdk-tests) sdk_tests="OFF"; shift ;;
        --jobs) jobs="${2:?--jobs needs a number}"; shift 2 ;;
        --build-dir) build_dir="${2:?--build-dir needs a path}"; shift 2 ;;
        --configure-only) configure_only=true; shift ;;
        --help|-h) usage ;;
        *) die "unknown argument: $1 (try --help)" ;;
    esac
done

# Resolve the prefix from the vendor variable when the caller did not name one.
if [ -z "$sdk_root" ]; then
    sdk_root="${ROCM_PATH:-${CUDAToolkit_ROOT:-${CUDA_HOME:-}}}"
fi

[ -n "$sdk_root" ] || die "no SDK root: set ONEFLOW_ACCEL_SDK_ROOT, export ROCM_PATH / CUDA_HOME, or pass --sdk-root <path>"
[ -d "$sdk_root" ] || die "SDK root is not a directory: $sdk_root"

# Infer the kind from the prefix, so the CMake call sets the matching enable flag.
if [ -x "$sdk_root/llvm/bin/clang++" ] || [ -d "$sdk_root/lib/cmake/hip" ]; then
    sdk="hip"
elif [ -x "$sdk_root/bin/nvcc" ]; then
    sdk="cuda"
else
    die "cannot identify the SDK under $sdk_root (no llvm/bin/clang++ or bin/nvcc)"
fi

case "$sdk" in
    hip)
        [ -x "$sdk_root/llvm/bin/clang++" ] || die "no HIP compiler under $sdk_root/llvm/bin/clang++"
        enable_var="-DONEFLOW_ENABLE_HIP=ON"
        tests_var="-DONEFLOW_ENABLE_HIP_TESTS=$sdk_tests"
        ;;
    cuda)
        [ -x "$sdk_root/bin/nvcc" ] || die "no nvcc under $sdk_root/bin/nvcc"
        enable_var="-DONEFLOW_ENABLE_CUDA=ON"
        tests_var="-DONEFLOW_ACCEL_BACKEND=CUDA"
        ;;
    *) die "unsupported SDK '$sdk' (expected hip or cuda)" ;;
esac

if [ -z "$jobs" ]; then
    jobs="$( (command -v nproc >/dev/null && nproc) || (command -v sysctl >/dev/null && sysctl -n hw.ncpu) || echo 4 )"
fi

cmake_args=(
    "$enable_var"
    "$tests_var"
    -DONEFLOW_ACCEL_SDK_ROOT="$sdk_root"
)

if [ -n "$arch" ]; then
    cmake_args+=( -DONEFLOW_ACCEL_SDK_ARCH="$arch" )
fi

if [ -n "$sdk_flags" ]; then
    cmake_args+=( -DONEFLOW_ACCEL_SDK_FLAGS="$sdk_flags" )
fi

printf 'SDK         : %s\n' "$sdk"
printf 'SDK root    : %s\n' "$sdk_root"
printf 'SDK arch    : %s\n' "${arch:-<node detection>}"
printf 'SDK flags   : %s\n' "${sdk_flags:-<none>}"
printf 'SDK tests   : %s\n' "$sdk_tests"
printf 'build dir   : %s\n' "$build_dir"
printf 'cmake       : %s\n' "$(command -v cmake)"
printf 'parallelism : %s\n' "$jobs"

cmake -S "$repo_root" -B "$build_dir" "${cmake_args[@]}"

if [ "$configure_only" = true ]; then
    exit 0
fi

cmake --build "$build_dir" --parallel "$jobs"
