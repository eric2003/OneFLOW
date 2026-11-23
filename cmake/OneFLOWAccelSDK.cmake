# Unified build-time accelerator SDK interface.
#
# Two independent axes:
#
#   ONEFLOW_ACCEL_BACKEND   runtime selection          (AUTO, CPU, HIP, CUDA, KOKKOS)
#   ONEFLOW_ACCEL_SDK_ROOT  build-time vendor SDK root (empty: configure no SDK)
#
# The SDK *kind* is inferred rather than asked for: from ONEFLOW_ENABLE_HIP /
# ONEFLOW_ENABLE_CUDA when they are set, otherwise from the contents of the root
# (llvm/bin/clang++ -> hip, bin/nvcc -> cuda). Only two optional knobs remain:
#
#   ONEFLOW_ACCEL_SDK_ROOT    prefix (DTK/ROCm for HIP, CUDA Toolkit for CUDA)
#   ONEFLOW_ACCEL_SDK_ARCH    gfx906 / sm_80      (empty: keep device detection)
#   ONEFLOW_ACCEL_SDK_FLAGS   extra SDK compiler flags (empty: none)
#
# Usage:
#
#     cmake -S . -B build -DONEFLOW_ENABLE_HIP=ON \
#           -DONEFLOW_ACCEL_SDK_ROOT=/path/to/dtk
#
#     cmake -S . -B build -DONEFLOW_ENABLE_CUDA=ON \
#           -DONEFLOW_ACCEL_SDK_ROOT=/path/to/cuda -DONEFLOW_ACCEL_SDK_ARCH=sm_80
#
# Nothing happens when the root is empty, so CPU and CI configurations are
# unaffected.
#
# Why a resolver is needed at all
# -------------------------------
# HIP/DCU SDKs do not configure themselves from a prefix:
#   * hip-config.cmake derives its import prefix from ROCM_PATH
#         set ( _IMPORT_PREFIX "${ROCM_PATH}/hip" "/opt/dtk/hip" )
#     so without ROCM_PATH it looks for libgalaxyhip under /opt/rocm and stops
#     with "the installation package was faulty ... but not all the files it
#     references";
#   * enable_language(HIP) runs find_package(amd_comgr) inside a try_compile
#     sub-project that only inherits the variables listed in
#     CMAKE_TRY_COMPILE_PLATFORM_VARIABLES.
# Both are derived here from ONEFLOW_ACCEL_SDK_ROOT.

set ( ONEFLOW_ACCEL_SDK_ROOT "" CACHE PATH
      "Accelerator SDK prefix (DTK/ROCm for hip, CUDA Toolkit for cuda). Empty: no SDK." )

set ( ONEFLOW_ACCEL_SDK_ARCH "" CACHE STRING
      "Target architecture(s) for the SDK, e.g. gfx906 or sm_80 (empty: keep detection)" )

set ( ONEFLOW_ACCEL_SDK_FLAGS "" CACHE STRING
      "Extra compiler flags for the SDK compiler, e.g. --gcc-toolchain=<host gcc> (empty: none)" )

function ( _OneFLOW_ResolveSDKKind out_var )
    if ( ONEFLOW_ENABLE_HIP OR ONEFLOW_HIP_EXPLICIT )
        set ( kind "hip" )
    elseif ( ONEFLOW_ENABLE_CUDA )
        set ( kind "cuda" )
    elseif ( EXISTS "${ONEFLOW_ACCEL_SDK_ROOT}/llvm/bin/clang++"
             OR IS_DIRECTORY "${ONEFLOW_ACCEL_SDK_ROOT}/lib/cmake/hip" )
        set ( kind "hip" )
    elseif ( EXISTS "${ONEFLOW_ACCEL_SDK_ROOT}/bin/nvcc" )
        set ( kind "cuda" )
    else ()
        message ( FATAL_ERROR
            "Cannot identify the accelerator SDK under '${ONEFLOW_ACCEL_SDK_ROOT}'. "
            "Set ONEFLOW_ENABLE_HIP/ONEFLOW_ENABLE_CUDA, or point "
            "ONEFLOW_ACCEL_SDK_ROOT at a DTK/ROCm or CUDA Toolkit prefix." )
    endif ()

    set ( ${out_var} "${kind}" PARENT_SCOPE )
endfunction ()

function ( _OneFLOW_ConfigureHipSDK root )
    message ( STATUS "OneFLOW accel SDK      = hip (${root})" )

    # hip-config.cmake takes its import prefix from here.
    set ( ENV{ROCM_PATH} "${root}" )
    set ( ROCM_PATH "${root}" CACHE PATH "ROCm/DTK prefix" FORCE )

    set ( sdk_prefixes )
    foreach ( candidate
              "${root}/dcc/lib64/cmake/AMDDeviceLibs"
              "${root}/dcc/comgr/lib64/cmake"
              "${root}"
              "${root}/.hyhal" )
        if ( IS_DIRECTORY "${candidate}" )
            list ( APPEND sdk_prefixes "${candidate}" )
        endif ()
    endforeach ()
    if ( sdk_prefixes )
        list ( REMOVE_DUPLICATES sdk_prefixes )
        set ( CMAKE_PREFIX_PATH "${sdk_prefixes};${CMAKE_PREFIX_PATH}"
              CACHE STRING "Package search prefixes" FORCE )
        message ( STATUS "OneFLOW accel prefixes = ${sdk_prefixes}" )
    endif ()

    # Dependencies hip-config.cmake requires, and the try_compile project needs
    # the same values, which is what CMAKE_TRY_COMPILE_PLATFORM_VARIABLES carries.
    set ( inherited "CMAKE_PREFIX_PATH" )
    foreach ( dep_pair
              "amd_comgr_DIR:${root}/dcc/comgr/lib64/cmake/amd_comgr"
              "AMDDeviceLibs_DIR:${root}/dcc/lib64/cmake/AMDDeviceLibs"
              "hsa-runtime64_DIR:${root}/.hyhal/lib/cmake/hsa-runtime64" )
        string ( REPLACE ":" ";" dep_pair_list "${dep_pair}" )
        list ( GET dep_pair_list 0 dep_name )
        list ( GET dep_pair_list 1 dep_dir )
        if ( IS_DIRECTORY "${dep_dir}" )
            set ( ${dep_name} "${dep_dir}" CACHE PATH
                  "Accelerator SDK dependency discovered by OneFLOWAccelSDK.cmake" FORCE )
            list ( APPEND inherited "${dep_name}" )
        endif ()
    endforeach ()

    set ( try_vars "${CMAKE_TRY_COMPILE_PLATFORM_VARIABLES};${inherited}" )
    list ( REMOVE_DUPLICATES try_vars )
    set ( CMAKE_TRY_COMPILE_PLATFORM_VARIABLES "${try_vars}" CACHE STRING
          "Variables forwarded into try_compile projects" FORCE )

    if ( NOT CMAKE_HIP_COMPILER AND EXISTS "${root}/llvm/bin/clang++" )
        set ( CMAKE_HIP_COMPILER "${root}/llvm/bin/clang++" CACHE FILEPATH
              "HIP compiler" FORCE )
        message ( STATUS "OneFLOW accel compiler = ${CMAKE_HIP_COMPILER}" )
    endif ()

    if ( NOT ONEFLOW_ACCEL_SDK_ARCH STREQUAL "" )
        set ( CMAKE_HIP_ARCHITECTURES "${ONEFLOW_ACCEL_SDK_ARCH}" CACHE STRING
              "HIP architectures" FORCE )
    endif ()

    if ( NOT ONEFLOW_ACCEL_SDK_FLAGS STREQUAL "" )
        set ( CMAKE_HIP_FLAGS "${ONEFLOW_ACCEL_SDK_FLAGS}" CACHE STRING
              "HIP compiler flags" FORCE )
    endif ()
endfunction ()

function ( _OneFLOW_ConfigureCudaSDK root )
    message ( STATUS "OneFLOW accel SDK      = cuda (${root})" )

    # find_package(CUDAToolkit) honours this instead of a PATH hunt.
    set ( CUDAToolkit_ROOT "${root}" CACHE PATH "CUDA Toolkit root" FORCE )

    if ( NOT CMAKE_CUDA_COMPILER AND EXISTS "${root}/bin/nvcc" )
        set ( CMAKE_CUDA_COMPILER "${root}/bin/nvcc" CACHE FILEPATH
              "CUDA compiler" FORCE )
        message ( STATUS "OneFLOW accel compiler = ${CMAKE_CUDA_COMPILER}" )
    endif ()

    if ( NOT ONEFLOW_ACCEL_SDK_ARCH STREQUAL "" )
        set ( CMAKE_CUDA_ARCHITECTURES "${ONEFLOW_ACCEL_SDK_ARCH}" CACHE STRING
              "CUDA architectures" FORCE )
    endif ()

    if ( NOT ONEFLOW_ACCEL_SDK_FLAGS STREQUAL "" )
        set ( CMAKE_CUDA_FLAGS "${ONEFLOW_ACCEL_SDK_FLAGS}" CACHE STRING
              "CUDA compiler flags" FORCE )
    endif ()
endfunction ()

function ( OneFLOW_ConfigureAccelSDK )
    set ( sdk_root "${ONEFLOW_ACCEL_SDK_ROOT}" )

    # Fall back to the vendor environment variable when a backend was requested
    # but no root was named.
    if ( sdk_root STREQUAL "" )
        if ( ( ONEFLOW_ENABLE_HIP OR ONEFLOW_HIP_EXPLICIT ) AND NOT "$ENV{ROCM_PATH}" STREQUAL "" )
            set ( sdk_root "$ENV{ROCM_PATH}" )
        elseif ( ONEFLOW_ENABLE_CUDA AND NOT "$ENV{CUDAToolkit_ROOT}" STREQUAL "" )
            set ( sdk_root "$ENV{CUDAToolkit_ROOT}" )
        elseif ( ONEFLOW_ENABLE_CUDA AND NOT "$ENV{CUDA_HOME}" STREQUAL "" )
            set ( sdk_root "$ENV{CUDA_HOME}" )
        endif ()
    endif ()

    if ( sdk_root STREQUAL "" )
        return()   # no SDK requested, nothing to do
    endif ()

    if ( NOT IS_DIRECTORY "${sdk_root}" )
        message ( FATAL_ERROR
            "ONEFLOW_ACCEL_SDK_ROOT is not a directory: ${sdk_root}" )
    endif ()

    _OneFLOW_ResolveSDKKind ( sdk_kind )
    if ( sdk_kind STREQUAL "hip" )
        _OneFLOW_ConfigureHipSDK ( "${sdk_root}" )
    else ()
        _OneFLOW_ConfigureCudaSDK ( "${sdk_root}" )
    endif ()
endfunction ()
