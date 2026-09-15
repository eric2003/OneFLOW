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
target_compile_features(solver_name_list_test PRIVATE cxx_std_17)
gtest_discover_tests(solver_name_list_test DISCOVERY_MODE PRE_TEST)


#set(CMAKE_FOLDER ${CMAKE_FOLDER}/solver)
#
#add_executable(solver_name_policy_test
#  solver_name_policy_test.cpp
#  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverMap.cpp
#  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverNameList.cpp
#  ${CMAKE_SOURCE_DIR}/codes/solver/src/Solver.cpp
#  ${CMAKE_SOURCE_DIR}/codes/solver/src/SolverState.cpp
#  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
#  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
#  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
#  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
#  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/CmxTask.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/Task.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/SimpleTask.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/Command.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/GridState.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/TaskState.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/ActionState.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/FileInfo.cpp
#  ${CMAKE_SOURCE_DIR}/codes/task/src/TaskRegister.cpp
#  ${CMAKE_SOURCE_DIR}/codes/global/src/Iteration.cpp
#  ${CMAKE_SOURCE_DIR}/codes/global/src/Ctrl.cpp
#  ${CMAKE_SOURCE_DIR}/codes/global/src/FileMap.cpp
#  ${CMAKE_SOURCE_DIR}/codes/zone/src/ZoneState.cpp
#  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
#  ${CMAKE_SOURCE_DIR}/codes/register/src/Register.cpp
#  ${CMAKE_SOURCE_DIR}/codes/register/src/Message.cpp
#  ${CMAKE_SOURCE_DIR}/codes/register/src/HXClone.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataPage.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataBook.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataPara.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataBase.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataField.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataBaseIO.cpp
#  ${CMAKE_SOURCE_DIR}/codes/database/src/DataBaseType.cpp
#  ${CMAKE_SOURCE_DIR}/codes/parallel/src/Parallel.cpp
#  ${CMAKE_SOURCE_DIR}/codes/parallel/src/BasicParallel.cpp
#)
#
#target_include_directories(solver_name_policy_test
#  PRIVATE
#    ${CMAKE_SOURCE_DIR}/codes/solver/include
#    ${CMAKE_SOURCE_DIR}/codes/basic/include
#	${CMAKE_SOURCE_DIR}/codes/task/include
#	${CMAKE_SOURCE_DIR}/codes/io/include
#	${CMAKE_SOURCE_DIR}/codes/global/include
#	${CMAKE_SOURCE_DIR}/codes/zone/include
#	${CMAKE_SOURCE_DIR}/codes/geometry/include
#	${CMAKE_SOURCE_DIR}/codes/register/include
#	${CMAKE_SOURCE_DIR}/codes/database/include
#	${CMAKE_SOURCE_DIR}/codes/array/include
#	${CMAKE_SOURCE_DIR}/codes/parallel/include
#)
#
#target_link_libraries(solver_name_policy_test PRIVATE GTest::gtest_main)
#target_compile_features(solver_name_policy_test PRIVATE cxx_std_17)
#
#gtest_discover_tests(solver_name_policy_test DISCOVERY_MODE PRE_TEST)