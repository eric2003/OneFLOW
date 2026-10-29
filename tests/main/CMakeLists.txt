set(CMAKE_FOLDER ${CMAKE_FOLDER}/main)

# ---------------------------------------------------------------------------
# TaskRegistry unit test
# ---------------------------------------------------------------------------
add_executable(task_registry_test
  task_registry_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuTask.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContext.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp
)

target_include_directories(task_registry_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/accel/include
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
# SimuContext + context-aware Task
# ---------------------------------------------------------------------------
add_executable(simu_context_test
  simu_context_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContext.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuTask.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateRegistry.cpp
  ${CMAKE_SOURCE_DIR}/codes/accel/src/EulerDomainStateLifecycle.cpp
)

target_include_directories(simu_context_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/accel/include
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


# ---------------------------------------------------------------------------
# Simulation command-line parser unit tests (independent of Prj and solver runtime)
# ---------------------------------------------------------------------------
add_executable(simulation_command_line_test
  simulation_command_line_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
)

target_include_directories(simulation_command_line_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(simulation_command_line_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(simulation_command_line_test PRIVATE cxx_std_20)

gtest_discover_tests(simulation_command_line_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# Prj compatibility adapter and path tests
# ---------------------------------------------------------------------------
add_executable(prj_cmdline_test
  prj_cmdline_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
)

target_include_directories(prj_cmdline_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(prj_cmdline_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(prj_cmdline_test PRIVATE cxx_std_20)

gtest_discover_tests(prj_cmdline_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# FieldPipeline stage metadata (header-only; no CFD / no FieldSimu.cpp)
# ---------------------------------------------------------------------------
add_executable(field_pipeline_test
  field_pipeline_test.cpp
)

target_include_directories(field_pipeline_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(field_pipeline_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(field_pipeline_test PRIVATE cxx_std_17)

gtest_discover_tests(field_pipeline_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# Optional control parameter lookup (isolated from full runtime initialization)
# ---------------------------------------------------------------------------
set(SIMU_CONTEXT_PARAMETER_DB_DIR ${CMAKE_SOURCE_DIR}/codes/database)
set(SIMU_CONTEXT_PARAMETER_BASIC_DIR ${CMAKE_SOURCE_DIR}/codes/basic)
set(SIMU_CONTEXT_PARAMETER_PARALLEL_DIR ${CMAKE_SOURCE_DIR}/codes/parallel)
set(SIMU_CONTEXT_PARAMETER_ZONE_DIR ${CMAKE_SOURCE_DIR}/codes/zone)

add_executable(simu_context_parameter_test
  simu_context_parameter_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextIO.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataBase.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataPara.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataBook.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataField.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataBaseType.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataPage.cpp
  ${SIMU_CONTEXT_PARAMETER_DB_DIR}/src/DataBaseIO.cpp
  ${SIMU_CONTEXT_PARAMETER_BASIC_DIR}/src/Fatal.cpp
  ${SIMU_CONTEXT_PARAMETER_PARALLEL_DIR}/src/BasicParallel.cpp
  ${SIMU_CONTEXT_PARAMETER_PARALLEL_DIR}/src/Parallel.cpp
  ${SIMU_CONTEXT_PARAMETER_ZONE_DIR}/src/ZoneState.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
)

target_include_directories(simu_context_parameter_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${SIMU_CONTEXT_PARAMETER_DB_DIR}/include
    ${SIMU_CONTEXT_PARAMETER_BASIC_DIR}/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
    ${SIMU_CONTEXT_PARAMETER_PARALLEL_DIR}/include
    ${SIMU_CONTEXT_PARAMETER_ZONE_DIR}/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
)

target_link_libraries(simu_context_parameter_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(simu_context_parameter_test PRIVATE cxx_std_20)

gtest_discover_tests(simu_context_parameter_test
  DISCOVERY_MODE PRE_TEST
)

# ---------------------------------------------------------------------------
# Case execution cleanup contract (header-only lifecycle boundary)
# ---------------------------------------------------------------------------
add_executable(simu_case_lifecycle_test
  simu_case_lifecycle_test.cpp
)

target_include_directories(simu_case_lifecycle_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/main/include
)

target_link_libraries(simu_case_lifecycle_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(simu_case_lifecycle_test PRIVATE cxx_std_17)

gtest_discover_tests(simu_case_lifecycle_test
  DISCOVERY_MODE PRE_TEST
)

