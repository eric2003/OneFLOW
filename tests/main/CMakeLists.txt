set(CMAKE_FOLDER ${CMAKE_FOLDER}/main)

# ---------------------------------------------------------------------------
# TaskRegistry unit test (no CFD kernel dependencies)
# ---------------------------------------------------------------------------
add_executable(task_registry_test
  task_registry_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuTask.cpp
)

target_include_directories(task_registry_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(task_registry_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(task_registry_test PRIVATE cxx_std_17)

gtest_discover_tests(task_registry_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# SimuBase interface smoke test
# ---------------------------------------------------------------------------
add_executable(simu_base_interface_test
  simu_base_interface_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuBase.cpp
)

target_include_directories(simu_base_interface_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(simu_base_interface_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(simu_base_interface_test PRIVATE cxx_std_17)

gtest_discover_tests(simu_base_interface_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# SimuContext pure logic (no MPI / control-file / accel)
# Links SimuContext.cpp + SimuTask.cpp only (not SimuContextEnv.cpp).
# ---------------------------------------------------------------------------
add_executable(simu_context_test
  simu_context_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContext.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuTask.cpp
)

target_include_directories(simu_context_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(simu_context_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(simu_context_test PRIVATE cxx_std_17)

gtest_discover_tests(simu_context_test
  DISCOVERY_MODE PRE_TEST
)
