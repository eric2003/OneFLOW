set(CMAKE_FOLDER ${CMAKE_FOLDER}/euler)

set(ONEFLOW_EULER_BACKEND_DIR
  ${CMAKE_SOURCE_DIR}/ports/kunshan/oneflow_1d_hip
)

add_executable(oneflow_euler_backend_cpu_test
  euler_backend_contract_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDWeno5.cpp
)

target_include_directories(oneflow_euler_backend_cpu_test
  PRIVATE
    ${ONEFLOW_EULER_BACKEND_DIR}
)

target_compile_options(oneflow_euler_backend_cpu_test
  PRIVATE
    -ffp-contract=off
)

target_link_libraries(oneflow_euler_backend_cpu_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(oneflow_euler_backend_cpu_test
  DISCOVERY_MODE PRE_TEST
)

add_executable(oneflow_flux_backend_bridge_test
  flux_backend_bridge_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDWeno5.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuFluxBackend.cpp
)

target_include_directories(oneflow_flux_backend_bridge_test
  PRIVATE
    ${ONEFLOW_EULER_BACKEND_DIR}
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CMAKE_SOURCE_DIR}/codes/project/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_compile_options(oneflow_flux_backend_bridge_test
  PRIVATE
    -ffp-contract=off
)

target_link_libraries(oneflow_flux_backend_bridge_test
  PRIVATE
)

# Bridge test is a standalone main(), not a gtest
# Run manually: ./bin/oneflow_flux_backend_bridge_test


add_executable(oneflow_euler_domain_contract_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_contract_test.cpp)
target_include_directories(oneflow_euler_domain_contract_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_contract_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_contract_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_state_registry_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_state_registry_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp)
target_include_directories(oneflow_euler_domain_state_registry_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_state_registry_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_state_registry_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_state_lifecycle_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_state_lifecycle_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp)
target_include_directories(oneflow_euler_domain_state_lifecycle_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_state_lifecycle_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_state_lifecycle_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_cpu_domain_backend_test
  ${CMAKE_SOURCE_DIR}/tests/euler_cpu_domain_backend_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuEulerDomainBackend.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp)
target_include_directories(oneflow_euler_cpu_domain_backend_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_cpu_domain_backend_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_cpu_domain_backend_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_domain_mrfield_adapter_test
  ${CMAKE_SOURCE_DIR}/tests/euler_domain_mrfield_adapter_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/EulerDomainMrFieldAdapter.cpp)
target_include_directories(oneflow_euler_domain_mrfield_adapter_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/main/include
  ${CMAKE_SOURCE_DIR}/codes/array/include
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_domain_mrfield_adapter_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_domain_mrfield_adapter_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_cpu_adapter_test
  ${CMAKE_SOURCE_DIR}/tests/euler_cpu_adapter_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerCpuAdapter.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/CpuFluxBackend.cpp)
target_include_directories(oneflow_euler_cpu_adapter_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_cpu_adapter_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_cpu_adapter_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_inv_flux_capability_test
  ${CMAKE_SOURCE_DIR}/tests/euler_inv_flux_capability_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerInvFluxCapability.cpp)
target_include_directories(oneflow_euler_inv_flux_capability_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_inv_flux_capability_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_inv_flux_capability_test DISCOVERY_MODE PRE_TEST)

add_executable(oneflow_euler_runge_kutta_capability_test
  ${CMAKE_SOURCE_DIR}/tests/euler_runge_kutta_capability_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/EulerRungeKuttaCapability.cpp)
target_include_directories(oneflow_euler_runge_kutta_capability_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/solver/include
  ${CMAKE_SOURCE_DIR}/codes/global/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
  ${CMAKE_SOURCE_DIR}/codes/project/include)
target_link_libraries(oneflow_euler_runge_kutta_capability_test PRIVATE GTest::gtest_main)
gtest_discover_tests(oneflow_euler_runge_kutta_capability_test DISCOVERY_MODE PRE_TEST)

# Hardware-dependent HIP/DCU contract. This is opt-in so ordinary CPU builds
# keep the root CTest suite reproducible on machines without a DCU. On a target
# node, configure with CMAKE_HIP_COMPILER and CMAKE_HIP_ARCHITECTURES, then the
# same target is discovered by the root CTest tree with the `HIP.` prefix and
# `hip`/`dcu` labels.
if (ONEFLOW_ENABLE_HIP_TESTS)
  if (NOT CMAKE_HIP_ARCHITECTURES)
    if (ONEFLOW_HIP_ARCHITECTURES)
      set(CMAKE_HIP_ARCHITECTURES "${ONEFLOW_HIP_ARCHITECTURES}"
          CACHE STRING "HIP/DCU architecture list" FORCE)
    else()
      message(FATAL_ERROR
        "ONEFLOW_ENABLE_HIP_TESTS=ON requires CMAKE_HIP_ARCHITECTURES "
        "or ONEFLOW_HIP_ARCHITECTURES")
    endif()
  endif()

  enable_language(HIP)
  find_package(hip CONFIG REQUIRED)
  include(${CMAKE_SOURCE_DIR}/cmake/OneFLOWEulerContract.cmake)
  oneflow_add_euler_hip_contract_test(
    oneflow_euler_backend_hip_test
    PORT_DIR ${ONEFLOW_EULER_BACKEND_DIR}
    ACCEL_DIR ${CMAKE_SOURCE_DIR}/codes/accel
    PROJECT_INC ${CMAKE_SOURCE_DIR}/codes/project/include
    BASIC_INC ${CMAKE_SOURCE_DIR}/codes/basic/include
    TEST_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/euler_backend_contract_test.cpp)

  add_test(
    NAME HIP.MainSolverFluxOneCall
    COMMAND $<TARGET_FILE:OneFLOWHipSmoke>)
  set_tests_properties(
    HIP.MainSolverFluxOneCall
    PROPERTIES LABELS "hardware;hip;dcu")
endif()
