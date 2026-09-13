cmake_minimum_required( VERSION 3.22 )
project( oneflow_1d_upwind LANGUAGES CXX )

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

include( ${CMAKE_CURRENT_LIST_DIR}/../../../cmake/UtilFunction.cmake )

add_executable( oneflow_1d_cpu CpuMain.cpp OneDUpwind.cpp )
add_executable(
    oneflow_1d_euler_cpu
    EulerMain.cpp
    OneDEuler.cpp
    OneDEulerBackend.cpp )
target_compile_options( oneflow_1d_euler_cpu PRIVATE -ffp-contract=off )
add_executable( oneflow_1d_weno5_dump Weno5DumpMain.cpp OneDWeno5.cpp )
target_compile_options( oneflow_1d_weno5_dump PRIVATE -ffp-contract=off )

option( ONEFLOW_1D_ENABLE_MPI_BENCHMARK "Build MPI Euler benchmarks" OFF )
if ( ONEFLOW_1D_ENABLE_MPI_BENCHMARK )
    find_package( MPI REQUIRED COMPONENTS CXX )
    add_executable(
        oneflow_1d_euler_mpi_cpu_benchmark
        EulerMpiBenchmark.cpp
        OneDEulerMpi.cpp )
    target_compile_options(
        oneflow_1d_euler_mpi_cpu_benchmark PRIVATE -ffp-contract=off )
    target_link_libraries(
        oneflow_1d_euler_mpi_cpu_benchmark PRIVATE MPI::MPI_CXX )
endif()


option(
    ONEFLOW_1D_ENABLE_GTEST
    "Build Euler backend GoogleTest targets"
    OFF )

if ( ONEFLOW_1D_ENABLE_GTEST )
    enable_testing()
    if ( NOT TARGET GTest::gtest_main )
        add_subdirectory(
            ${CMAKE_CURRENT_LIST_DIR}/../../../ThirdParty/googletest
            ${CMAKE_CURRENT_BINARY_DIR}/googletest
            EXCLUDE_FROM_ALL )
    endif()
    include( GoogleTest )

    set( ONEFLOW_EULER_CONTRACT_TEST
         ${CMAKE_CURRENT_LIST_DIR}/../../../tests/euler/euler_backend_contract_test.cpp )
    add_executable(
        oneflow_1d_euler_cpu_backend_test
        ${ONEFLOW_EULER_CONTRACT_TEST}
        OneDEuler.cpp
        OneDEulerBackend.cpp
        OneDWeno5.cpp )
    target_include_directories(
        oneflow_1d_euler_cpu_backend_test
        PRIVATE ${CMAKE_CURRENT_LIST_DIR} )
    target_compile_options(
        oneflow_1d_euler_cpu_backend_test PRIVATE -ffp-contract=off )
    target_link_libraries(
        oneflow_1d_euler_cpu_backend_test PRIVATE GTest::gtest_main )
    gtest_discover_tests(
            oneflow_1d_euler_cpu_backend_test
            DISCOVERY_MODE PRE_TEST
            TEST_PREFIX "CPU." )
endif()


# ---------------------------------------------------------------------------
# OneFLOW accelerator backend infrastructure (codes/accel)
# Shared between port and main solver; provides AccelRuntime, AccelBackend,
# DeviceBuffer<T>, and the HIP/CUDA/Kokkos backend registry.
# ---------------------------------------------------------------------------
set( ONEFLOW_ACCEL_DIR ${CMAKE_CURRENT_LIST_DIR}/../../../codes/accel )
set( ONEFLOW_PROJECT_INC ${CMAKE_CURRENT_LIST_DIR}/../../../codes/project/include )
set( ONEFLOW_BASIC_INC   ${CMAKE_CURRENT_LIST_DIR}/../../../codes/basic/include )

set( ONEFLOW_ACCEL_SOURCES
    ${ONEFLOW_ACCEL_DIR}/src/AccelBackend.cpp
    ${ONEFLOW_ACCEL_DIR}/src/AccelRuntime.cpp
    ${ONEFLOW_ACCEL_DIR}/src/CpuBackend.cpp
)

set( ONEFLOW_ACCEL_HIP_SOURCES
    ${ONEFLOW_ACCEL_DIR}/src/HipBackend.cpp
    ${ONEFLOW_ACCEL_DIR}/src/HipKernel.hip
)

set( ONEFLOW_ACCEL_INCLUDES
    ${ONEFLOW_ACCEL_DIR}/include
    ${ONEFLOW_PROJECT_INC}
    ${ONEFLOW_BASIC_INC}
)

set( ONEFLOW_ACCEL_DEFINES
    ONEFLOW_ENABLE_HIP
    ONEFLOW_ENABLE_MULTI_DEVICE
    ONEFLOW_DEFAULT_ACCEL_BACKEND="HIP"
)

# Helper: attach accel backend to a HIP target
macro( oneflow_accel_attach target )
    target_sources( ${target} PRIVATE
        ${ONEFLOW_ACCEL_SOURCES}
        ${ONEFLOW_ACCEL_HIP_SOURCES}
    )
    target_include_directories( ${target} PRIVATE ${ONEFLOW_ACCEL_INCLUDES} )
    target_compile_definitions( ${target} PRIVATE ${ONEFLOW_ACCEL_DEFINES} )
    target_link_libraries( ${target} PRIVATE hip::host )
endmacro()


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
    oneflow_accel_attach( oneflow_1d_hip )
    set_source_files_properties( OneDEuler.hip PROPERTIES LANGUAGE HIP )
    set_source_files_properties( OneDEulerPersistent.hip PROPERTIES LANGUAGE HIP )
    add_executable(
        oneflow_1d_euler_hip
        EulerMain.cpp
        OneDEuler.cpp
        OneDEuler.hip
        OneDEulerBackend.cpp
        OneDEulerPersistent.hip )
    target_compile_definitions( oneflow_1d_euler_hip PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_euler_hip PRIVATE -ffp-contract=off )
    oneflow_accel_attach( oneflow_1d_euler_hip )

    if ( ONEFLOW_1D_ENABLE_GTEST )
        add_executable(
            oneflow_1d_euler_hip_backend_test
            ${ONEFLOW_EULER_CONTRACT_TEST}
            OneDEuler.cpp
            OneDEulerBackend.cpp
            OneDEulerPersistent.hip )
        target_compile_definitions(
            oneflow_1d_euler_hip_backend_test PRIVATE ONEFLOW_1D_USE_HIP=1 )
        target_compile_options(
            oneflow_1d_euler_hip_backend_test PRIVATE -ffp-contract=off )
        target_include_directories(
            oneflow_1d_euler_hip_backend_test
            PRIVATE ${CMAKE_CURRENT_LIST_DIR} )
        oneflow_accel_attach( oneflow_1d_euler_hip_backend_test )
        target_link_libraries(
            oneflow_1d_euler_hip_backend_test
            PRIVATE GTest::gtest_main )
        gtest_discover_tests(
            oneflow_1d_euler_hip_backend_test
            DISCOVERY_MODE PRE_TEST
            TEST_PREFIX "HIP." )
    endif()


    set_source_files_properties( OneDEulerProfile.hip PROPERTIES LANGUAGE HIP )
    add_executable(
        oneflow_1d_euler_benchmark
        EulerBenchmark.cpp
        OneDEuler.cpp
        OneDEulerBackend.cpp
        OneDEulerProfile.hip )
    target_compile_options( oneflow_1d_euler_benchmark PRIVATE -ffp-contract=off )
    oneflow_accel_attach( oneflow_1d_euler_benchmark )

    add_executable( oneflow_1d_euler_stateful_benchmark EulerStatefulBenchmark.cpp OneDEuler.cpp OneDEulerBackend.cpp OneDEulerPersistent.hip )
    target_compile_definitions( oneflow_1d_euler_stateful_benchmark PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_euler_stateful_benchmark PRIVATE -ffp-contract=off )
    oneflow_accel_attach( oneflow_1d_euler_stateful_benchmark )
    add_executable( oneflow_1d_euler_gpu_only_benchmark EulerGpuOnlyBenchmark.cpp OneDEuler.cpp OneDEulerBackend.cpp OneDEulerPersistent.hip )
    target_compile_definitions( oneflow_1d_euler_gpu_only_benchmark PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_euler_gpu_only_benchmark PRIVATE -ffp-contract=off )
    oneflow_accel_attach( oneflow_1d_euler_gpu_only_benchmark )

    if ( ONEFLOW_1D_ENABLE_MPI_BENCHMARK )
        set_source_files_properties( OneDEulerMpi.hip PROPERTIES LANGUAGE HIP )
        add_executable(
            oneflow_1d_euler_mpi_hip_benchmark
            EulerMpiBenchmark.cpp
            OneDEulerMpi.cpp
            OneDEulerMpi.hip )
        target_compile_definitions(
            oneflow_1d_euler_mpi_hip_benchmark
            PRIVATE ONEFLOW_1D_MPI_USE_HIP=1 )
        target_compile_options(
            oneflow_1d_euler_mpi_hip_benchmark PRIVATE -ffp-contract=off )
        oneflow_accel_attach( oneflow_1d_euler_mpi_hip_benchmark )
        target_link_libraries(
            oneflow_1d_euler_mpi_hip_benchmark
            PRIVATE MPI::MPI_CXX )

        add_executable(
            oneflow_1d_euler_mpi_regression
            EulerMpiRegression.cpp
            OneDEulerMpi.cpp
            OneDEulerMpi.hip )
        target_compile_definitions(
            oneflow_1d_euler_mpi_regression
            PRIVATE ONEFLOW_1D_MPI_USE_HIP=1 )
        target_compile_options(
            oneflow_1d_euler_mpi_regression PRIVATE -ffp-contract=off )
        oneflow_accel_attach( oneflow_1d_euler_mpi_regression )
        target_link_libraries(
            oneflow_1d_euler_mpi_regression
            PRIVATE MPI::MPI_CXX )
    endif()

    set_source_files_properties( OneDWeno5.hip PROPERTIES LANGUAGE HIP )
    add_executable( oneflow_1d_weno5_replay Weno5ReplayMain.cpp OneDWeno5.cpp OneDWeno5.hip )
    target_compile_definitions( oneflow_1d_weno5_replay PRIVATE ONEFLOW_1D_USE_HIP=1 )
    target_compile_options( oneflow_1d_weno5_replay PRIVATE -ffp-contract=off )
    oneflow_accel_attach( oneflow_1d_weno5_replay )
endif()
