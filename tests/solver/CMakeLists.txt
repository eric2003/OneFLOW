set(CMAKE_FOLDER ${CMAKE_FOLDER}/solver)

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