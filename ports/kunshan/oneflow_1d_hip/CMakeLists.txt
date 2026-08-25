cmake_minimum_required( VERSION 3.22 )
project( oneflow_1d_upwind LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

add_executable( oneflow_1d_cpu CpuMain.cpp OneDUpwind.cpp )
add_executable( oneflow_1d_euler_cpu EulerMain.cpp OneDEuler.cpp )
target_compile_options( oneflow_1d_euler_cpu PRIVATE -ffp-contract=off )
add_executable( oneflow_1d_weno5_dump Weno5DumpMain.cpp OneDWeno5.cpp )
target_compile_options( oneflow_1d_weno5_dump PRIVATE -ffp-contract=off )

option( ONEFLOW_1D_ENABLE_HIP "Build the native HIP validation executable" OFF )
if ( ONEFLOW_1D_ENABLE_HIP )
    enable_language( HIP )
    find_package( hip CONFIG REQUIRED )
    if ( NOT DEFINED CMAKE_HIP_ARCHITECTURES )
        set( CMAKE_HIP_ARCHITECTURES gfx906 )
    endif()
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
