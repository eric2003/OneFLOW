cmake_minimum_required( VERSION 3.22 )
project( oneflow_1d_upwind LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

include( ${CMAKE_CURRENT_LIST_DIR}/../../../cmake/UtilFunction.cmake )

add_executable( oneflow_1d_cpu CpuMain.cpp OneDUpwind.cpp )
add_executable( oneflow_1d_euler_cpu EulerMain.cpp OneDEuler.cpp )
target_compile_options( oneflow_1d_euler_cpu PRIVATE -ffp-contract=off )
add_executable( oneflow_1d_weno5_dump Weno5DumpMain.cpp OneDWeno5.cpp )
target_compile_options( oneflow_1d_weno5_dump PRIVATE -ffp-contract=off )

option( ONEFLOW_1D_ENABLE_HIP "Build the native HIP validation executable" OFF )
if ( ONEFLOW_1D_ENABLE_HIP )
    if ( NOT CMAKE_HIP_ARCHITECTURES )
        DetectHIPArchitectures( ONEFLOW_1D_HIP_ARCHITECTURES )
        if ( ONEFLOW_1D_HIP_ARCHITECTURES )
            set( CMAKE_HIP_ARCHITECTURES "${ONEFLOW_1D_HIP_ARCHITECTURES}"
                 CACHE STRING "HIP/DCU architecture list" FORCE )
        else()
            message( FATAL_ERROR
                "ONEFLOW_1D_ENABLE_HIP=ON but no HIP architecture was detected. "
                "Set CMAKE_HIP_ARCHITECTURES or ONEFLOW_HIP_ARCHITECTURES "
                "for cross-compilation." )
        endif()
    endif()
    enable_language( HIP )
    find_package( hip CONFIG REQUIRED )
    set_source_files_properties( OneDUpwind.hip PROPERTIES LANGUAGE HIP )
    add_executable(
        oneflow_1d_hip
        HipMain.cpp
        OneDUpwind.cpp
        OneDUpwind.hip )
    target_link_libraries( oneflow_1d_hip PRIVATE hip::host )
    set_source_files_properties( OneDEuler.hip PROPERTIES LANGUAGE HIP )
    add_executable( oneflow_1d_euler_hip EulerMain.cpp OneDEuler.cpp OneDEuler.hip )
    target_compile_definitions( oneflow_1d_euler_hip PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_euler_hip PRIVATE -ffp-contract=off )
    target_link_libraries( oneflow_1d_euler_hip PRIVATE hip::host )

    set_source_files_properties( OneDWeno5.hip PROPERTIES LANGUAGE HIP )
    add_executable( oneflow_1d_weno5_replay Weno5ReplayMain.cpp OneDWeno5.cpp OneDWeno5.hip )
    target_compile_definitions( oneflow_1d_weno5_replay PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_weno5_replay PRIVATE -ffp-contract=off )
    target_link_libraries( oneflow_1d_weno5_replay PRIVATE hip::host )
endif()
