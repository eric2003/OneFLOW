set(CMAKE_FOLDER ${CMAKE_FOLDER}/accel)

# Standalone accounting contracts do not require a solver or accelerator.
add_executable(stage_profiler_test
  stage_profiler_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/StageProfiler.cpp
)
target_include_directories(stage_profiler_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/accel/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
)
target_link_libraries(stage_profiler_test PRIVATE GTest::gtest_main)
target_compile_features(stage_profiler_test PRIVATE cxx_std_17)
gtest_discover_tests(stage_profiler_test DISCOVERY_MODE PRE_TEST)
add_test(NAME StageProfiler.HostAccountingExplicitlyDisabled
  COMMAND stage_profiler_test --gtest_filter=StageProfiler.SwitchControlsSidecar)
set_tests_properties(StageProfiler.HostAccountingExplicitlyDisabled PROPERTIES
  ENVIRONMENT "ONEFLOW_HOST_ACCOUNTING=0")

add_executable(accel_runtime_test
  accel_runtime_test.cpp
)

target_link_libraries(accel_runtime_test
  PRIVATE
    OneFLOW_RuntimeObjects
    GTest::gtest_main
)

target_compile_features(accel_runtime_test PRIVATE cxx_std_17)

gtest_discover_tests(accel_runtime_test
  DISCOVERY_MODE PRE_TEST
)
