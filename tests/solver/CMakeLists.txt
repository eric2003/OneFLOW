set(CMAKE_FOLDER ${CMAKE_FOLDER}/solver)

# Pure policy only - no SolverMap, no IO, no parallel.
add_executable(solver_name_policy_test
  solver_name_policy_test.cpp
)

target_include_directories(solver_name_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(solver_name_policy_test PRIVATE GTest::gtest_main)
target_compile_features(solver_name_policy_test PRIVATE cxx_std_17)
gtest_discover_tests(solver_name_policy_test DISCOVERY_MODE PRE_TEST)


# Name list + policy + minimal IO (file parse path optional).
add_executable(solver_name_list_test
  solver_name_list_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverNameList.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
)

target_include_directories(solver_name_list_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/io/include
	${CMAKE_SOURCE_DIR}/codes/task/include
)

target_link_libraries(solver_name_list_test PRIVATE GTest::gtest_main)
target_compile_features(solver_name_list_test PRIVATE cxx_std_20)
gtest_discover_tests(solver_name_list_test DISCOVERY_MODE PRE_TEST)


# Index maps + SelectSolverNames - no SafeClone / SolverState / Lusgs.
add_executable(solver_map_index_test
  solver_map_index_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverMapIndex.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverNameList.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
)

target_include_directories(solver_map_index_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/task/include
)

target_link_libraries(solver_map_index_test PRIVATE GTest::gtest_main)
target_compile_features(solver_map_index_test PRIVATE cxx_std_20)
gtest_discover_tests(solver_map_index_test DISCOVERY_MODE PRE_TEST)



# Pure field-name contracts (no Zone / SolverState / DataBase).
add_executable(field_name_contract_test
  field_name_contract_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/FieldNameContract.cpp
)

target_include_directories(field_name_contract_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
)

target_link_libraries(field_name_contract_test PRIVATE GTest::gtest_main)
target_compile_features(field_name_contract_test PRIVATE cxx_std_17)
gtest_discover_tests(field_name_contract_test DISCOVERY_MODE PRE_TEST)


# Per-solver field-name tables (ns/ins/turb modules; still no Zone/DataBase).
add_executable(solver_field_names_test
  solver_field_names_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/FieldNameContract.cpp
  ${CMAKE_SOURCE_DIR}/codes/ns/src/NsFieldNames.cpp
  ${CMAKE_SOURCE_DIR}/codes/ins/src/INsFieldNames.cpp
  ${CMAKE_SOURCE_DIR}/codes/turb/src/TurbFieldNames.cpp
)

target_include_directories(solver_field_names_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/ns/include
    ${CMAKE_SOURCE_DIR}/codes/ins/include
    ${CMAKE_SOURCE_DIR}/codes/turb/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
)

target_link_libraries(solver_field_names_test PRIVATE GTest::gtest_main)
target_compile_features(solver_field_names_test PRIVATE cxx_std_17)
gtest_discover_tests(solver_field_names_test DISCOVERY_MODE PRE_TEST)


# Mean-flow cross-solver view contract (no Zone/DataBase bootstrap).
add_executable(mean_flow_view_test
  mean_flow_view_test.cpp
)

target_include_directories(mean_flow_view_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
)

target_link_libraries(mean_flow_view_test PRIVATE GTest::gtest_main)
target_compile_features(mean_flow_view_test PRIVATE cxx_std_17)
gtest_discover_tests(mean_flow_view_test DISCOVERY_MODE PRE_TEST)


# Zone-level shared field name contract (no Zone bootstrap required for names).
add_executable(zone_shared_fields_test
  zone_shared_fields_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/ZoneSharedFields.cpp
)

target_include_directories(zone_shared_fields_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/solver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
)

target_link_libraries(zone_shared_fields_test PRIVATE GTest::gtest_main)
target_compile_features(zone_shared_fields_test PRIVATE cxx_std_17)
gtest_discover_tests(zone_shared_fields_test DISCOVERY_MODE PRE_TEST)


# Pure scalar limiter policies; no mesh, solver state, or runtime bootstrap.
add_executable(limiter_policy_test
  limiter_policy_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/usolver/src/ULimiterPolicy.cpp
)

target_include_directories(limiter_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/usolver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(limiter_policy_test PRIVATE GTest::gtest_main)
target_compile_features(limiter_policy_test PRIVATE cxx_std_20)
gtest_discover_tests(limiter_policy_test DISCOVERY_MODE PRE_TEST)

 
# Neighbor-range construction tests use synthetic face connectivity without solver startup.
add_executable(limiter_range_test
  limiter_range_test.cpp
)

target_include_directories(limiter_range_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/usolver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(limiter_range_test PRIVATE GTest::gtest_main)
target_compile_features(limiter_range_test PRIVATE cxx_std_20)
gtest_discover_tests(limiter_range_test DISCOVERY_MODE PRE_TEST)


# Pure gradient face policies; no mesh or solver startup.
add_executable(grad_policy_test
  grad_policy_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/usolver/src/UGradPolicy.cpp
)

target_include_directories(grad_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/usolver/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(grad_policy_test PRIVATE GTest::gtest_main)
target_compile_features(grad_policy_test PRIVATE cxx_std_20)
gtest_discover_tests(grad_policy_test DISCOVERY_MODE PRE_TEST)

# Pure Euler convective flux helpers; no mesh or solver bootstrap.
add_executable(inv_flux_policy_test
  inv_flux_policy_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/ns/src/NsInvFluxPolicy.cpp
)

target_include_directories(inv_flux_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/ns/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(inv_flux_policy_test PRIVATE GTest::gtest_main)
target_compile_features(inv_flux_policy_test PRIVATE cxx_std_20)
gtest_discover_tests(inv_flux_policy_test DISCOVERY_MODE PRE_TEST)

add_executable(vis_flux_policy_test
  vis_flux_policy_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/ns/src/NsVisFluxPolicy.cpp
)
target_include_directories(vis_flux_policy_test PRIVATE
  ${CMAKE_SOURCE_DIR}/codes/ns/include
  ${CMAKE_SOURCE_DIR}/codes/basic/include
)
target_link_libraries(vis_flux_policy_test PRIVATE GTest::gtest_main)
target_compile_features(vis_flux_policy_test PRIVATE cxx_std_20)
gtest_discover_tests(vis_flux_policy_test DISCOVERY_MODE PRE_TEST)

# Pure Sutherland viscosity formulas; no mesh bootstrap.
add_executable(sutherland_policy_test
  sutherland_policy_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/model/src/Sutherland.cpp
)

target_include_directories(sutherland_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/model/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(sutherland_policy_test PRIVATE GTest::gtest_main)
target_compile_features(sutherland_policy_test PRIVATE cxx_std_20)
gtest_discover_tests(sutherland_policy_test DISCOVERY_MODE PRE_TEST)
