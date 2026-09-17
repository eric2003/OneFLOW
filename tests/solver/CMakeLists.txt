set(CMAKE_FOLDER ${CMAKE_FOLDER}/solver)

# Pure policy only — no SolverMap, no IO, no parallel.
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


# Index maps + SelectSolverNames — no SafeClone / SolverState / Lusgs.
add_executable(solver_map_index_test
  solver_map_index_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverMapIndex.cpp
  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverNameList.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
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

