set(CMAKE_FOLDER ${CMAKE_FOLDER}/task)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
set(TASK_INCLUDE_DIR ${OneFLOW_Root}/task/include)
set(TASK_SRC_DIR ${OneFLOW_Root}/task/src)
set(IO_INCLUDE_DIR ${OneFLOW_Root}/io/include)
set(IO_SRC_DIR ${OneFLOW_Root}/io/src)
set(SOLVER_INCLUDE_DIR ${OneFLOW_Root}/solver/include)
set(SOLVER_SRC_DIR ${OneFLOW_Root}/solver/src)
set(ZONE_INCLUDE_DIR ${OneFLOW_Root}/zone/include)
set(ZONE_SRC_DIR ${OneFLOW_Root}/zone/src)
set(GEOMETRY_INCLUDE_DIR ${OneFLOW_Root}/geometry/include)
set(GLOBAL_INCLUDE_DIR ${OneFLOW_Root}/global/include)
set(GLOBAL_SRC_DIR ${OneFLOW_Root}/global/src)
set(REGISTER_INCLUDE_DIR ${OneFLOW_Root}/register/include)
set(REGISTER_SRC_DIR ${OneFLOW_Root}/register/src)
set(DATABASE_INCLUDE_DIR ${OneFLOW_Root}/database/include)
set(DATABASE_SRC_DIR ${OneFLOW_Root}/database/src)
set(PARALLEL_INCLUDE_DIR ${OneFLOW_Root}/parallel/include)
set(PARALLEL_SRC_DIR ${OneFLOW_Root}/parallel/src)
set(ARRAY_INCLUDE_DIR ${OneFLOW_Root}/array/include)

set(DATABASE_TEST_SOURCES
    action_map_test.cpp
	ActionMapReadFileTest.cpp
	TaskRegisterTest.cpp
	category_test.cpp
	CmxTaskDispatchTest.cpp
	TaskTest.cpp
	command_test.cpp
	cmd_test.cpp
	${TASK_SRC_DIR}/ActionMap.cpp
	${TASK_SRC_DIR}/Category.cpp
	${TASK_SRC_DIR}/TaskRegister.cpp
	${TASK_SRC_DIR}/Command.cpp
	${TASK_SRC_DIR}/Task.cpp
	${TASK_SRC_DIR}/TaskState.cpp
	${TASK_SRC_DIR}/SimpleTask.cpp
	${TASK_SRC_DIR}/CmxTask.cpp	
	${TASK_SRC_DIR}/FileInfo.cpp
	${TASK_SRC_DIR}/ActionState.cpp
	${TASK_SRC_DIR}/GridState.cpp
	${REGISTER_SRC_DIR}/Message.cpp
	${REGISTER_SRC_DIR}/Register.cpp
	${REGISTER_SRC_DIR}/HXClone.cpp
	${IO_SRC_DIR}/TextFileParser.cpp
	${IO_SRC_DIR}/CommentLine.cpp
	${IO_SRC_DIR}/Prj.cpp
	${IO_SRC_DIR}/FileUtils.cpp
	${BASIC_SRC_DIR}/Word.cpp
	${SOLVER_SRC_DIR}/SolverState.cpp
	${SOLVER_SRC_DIR}/SolverMap.cpp
	${SOLVER_SRC_DIR}/Solver.cpp
	${DATABASE_SRC_DIR}/DataPage.cpp
	${DATABASE_SRC_DIR}/DataBook.cpp
	${DATABASE_SRC_DIR}/DataPara.cpp
	${DATABASE_SRC_DIR}/DataBaseType.cpp
	${DATABASE_SRC_DIR}/DataBase.cpp
	${DATABASE_SRC_DIR}/DataField.cpp
	${DATABASE_SRC_DIR}/DataBaseIO.cpp
	${ZONE_SRC_DIR}/ZoneState.cpp
	${GLOBAL_SRC_DIR}/Iteration.cpp
	${GLOBAL_SRC_DIR}/Ctrl.cpp
	${GLOBAL_SRC_DIR}/FileMap.cpp
	${PARALLEL_SRC_DIR}/BasicParallel.cpp
	${PARALLEL_SRC_DIR}/Parallel.cpp
	
	
	
)

add_executable(task_test
  ${DATABASE_TEST_SOURCES}
)

target_include_directories(task_test
  PRIVATE
    ${TASK_INCLUDE_DIR}
	${BASIC_INCLUDE_DIR}
	${IO_INCLUDE_DIR}
	${SOLVER_INCLUDE_DIR}
	${REGISTER_INCLUDE_DIR}
	${DATABASE_INCLUDE_DIR}
	${ZONE_INCLUDE_DIR}
	${GLOBAL_INCLUDE_DIR}
	${GEOMETRY_INCLUDE_DIR}
	${PARALLEL_INCLUDE_DIR}
	${ARRAY_INCLUDE_DIR}
	
)

target_link_libraries(task_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(task_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
