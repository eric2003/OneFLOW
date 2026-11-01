set(CMAKE_FOLDER ${CMAKE_FOLDER}/accel)

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
