include_guard( GLOBAL )

# Register the standalone/target-node Euler HIP contract in the same CTest
# topology as the root project. The caller owns HIP language initialization,
# architecture selection, and the GTest target.
function( oneflow_add_euler_hip_contract_test target )
    set( options )
    set( one_value_args PORT_DIR ACCEL_DIR PROJECT_INC BASIC_INC TEST_SOURCE )
    cmake_parse_arguments(
        ONEFLOW_EULER_HIP
        "${options}"
        "${one_value_args}"
        ""
        ${ARGN} )

    foreach( required PORT_DIR ACCEL_DIR PROJECT_INC BASIC_INC TEST_SOURCE )
        if( NOT ONEFLOW_EULER_HIP_${required} )
            message( FATAL_ERROR
                "oneflow_add_euler_hip_contract_test requires ${required}" )
        endif()
    endforeach()

    set_source_files_properties(
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDWeno5.hip
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDEulerPersistent.hip
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/HipKernel.hip
        PROPERTIES LANGUAGE HIP )

    add_executable(
        ${target}
        ${ONEFLOW_EULER_HIP_TEST_SOURCE}
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDEuler.cpp
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDEulerBackend.cpp
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDWeno5.cpp
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDWeno5.hip
        ${ONEFLOW_EULER_HIP_PORT_DIR}/OneDEulerPersistent.hip
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/AccelBackend.cpp
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/AccelRuntime.cpp
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/CpuBackend.cpp
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/HipBackend.cpp
        ${ONEFLOW_EULER_HIP_ACCEL_DIR}/src/HipKernel.hip )
    target_compile_definitions(
        ${target}
        PRIVATE
            ONEFLOW_1D_USE_HIP=1
            ONEFLOW_ENABLE_HIP
            ONEFLOW_ENABLE_MULTI_DEVICE
            ONEFLOW_DEFAULT_ACCEL_BACKEND="HIP" )
    target_compile_options( ${target} PRIVATE -ffp-contract=off )
    target_include_directories(
        ${target}
        PRIVATE
            ${ONEFLOW_EULER_HIP_PORT_DIR}
            ${ONEFLOW_EULER_HIP_ACCEL_DIR}/include
            ${ONEFLOW_EULER_HIP_PROJECT_INC}
            ${ONEFLOW_EULER_HIP_BASIC_INC} )
    target_link_libraries( ${target} PRIVATE GTest::gtest_main hip::host )
    gtest_discover_tests(
        ${target}
        DISCOVERY_MODE PRE_TEST
        TEST_PREFIX "HIP."
        PROPERTIES LABELS "hardware;hip;dcu" )
endfunction()
