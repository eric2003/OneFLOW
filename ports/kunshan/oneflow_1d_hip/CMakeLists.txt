cmake_minimum_required( VERSION 3.22 )
project( oneflow_1d_upwind LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

add_executable( oneflow_1d_cpu CpuMain.cpp OneDUpwind.cpp )

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
endif()
