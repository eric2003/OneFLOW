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

# Hardware-dependent HIP/DCU contract. This is opt-in so ordinary CPU builds
# keep the root CTest suite reproducible on machines without a DCU. On a target
# node, configure with CMAKE_HIP_COMPILER and CMAKE_HIP_ARCHITECTURES, then the
# same target is discovered by the root CTest tree with the `HIP.` prefix and
# `hardware`/`hip`/`dcu` labels.
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
endif()
