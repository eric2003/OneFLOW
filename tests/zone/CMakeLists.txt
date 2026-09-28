set(CMAKE_FOLDER ${CMAKE_FOLDER}/zone)

add_executable(wall_dist_policy_test
  wall_dist_policy_test.cpp
)

target_include_directories(wall_dist_policy_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/zone/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(wall_dist_policy_test PRIVATE GTest::gtest_main)
target_compile_features(wall_dist_policy_test PRIVATE cxx_std_20)

gtest_discover_tests(wall_dist_policy_test DISCOVERY_MODE PRE_TEST)