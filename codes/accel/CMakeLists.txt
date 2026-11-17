if ( NOT ONEFLOW_BUILD_HIP )
    return()
endif()

add_executable ( OneFLOWHipSmoke
    src/HipSmoke.cpp
    src/HipKernel.hip
    src/AccelBackend.cpp
    src/AccelRuntime.cpp
    src/CpuBackend.cpp
    src/CpuFluxBackend.cpp
    src/EulerCpuAdapter.cpp
    src/HipBackend.cpp
    src/HipFluxBackend.hip
    src/StageProfiler.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BoundaryExchangeContract.cpp
)

target_include_directories ( OneFLOWHipSmoke
    PRIVATE
        include
        ${CMAKE_SOURCE_DIR}/codes/project/include
        ${CMAKE_SOURCE_DIR}/codes/basic/include
        ${CMAKE_SOURCE_DIR}/codes/parallel/include
        ${hip_INCLUDE_DIRS}
)

target_compile_definitions ( OneFLOWHipSmoke
    PRIVATE
        ONEFLOW_ENABLE_HIP
        ONEFLOW_ENABLE_MULTI_DEVICE
        ONEFLOW_DEFAULT_ACCEL_BACKEND="HIP"
)

target_link_libraries ( OneFLOWHipSmoke PRIVATE hip::host )
