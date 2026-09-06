
set(ONEFLOW_EULER_BACKEND_DIR
  ${CMAKE_SOURCE_DIR}/ports/kunshan/oneflow_1d_hip
)

add_executable(oneflow_euler_backend_cpu_test
  euler_backend_contract_test.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEuler.cpp
  ${ONEFLOW_EULER_BACKEND_DIR}/OneDEulerBackend.cpp
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
