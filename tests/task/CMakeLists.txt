set(CMAKE_FOLDER ${CMAKE_FOLDER}/task)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
set(TASK_INCLUDE_DIR ${OneFLOW_Root}/task/include)
set(TASK_SRC_DIR ${OneFLOW_Root}/task/src)
set(IO_INCLUDE_DIR ${OneFLOW_Root}/io/include)
set(IO_SRC_DIR ${OneFLOW_Root}/io/src)


set(DATABASE_TEST_SOURCES
    action_map_test.cpp
	ActionMapReadFileTest.cpp
	TaskRegisterTest.cpp
	category_test.cpp
	${TASK_SRC_DIR}/ActionMap.cpp
	${TASK_SRC_DIR}/Category.cpp
	${TASK_SRC_DIR}/TaskRegister.cpp
	${IO_SRC_DIR}/TextFileParser.cpp
	${IO_SRC_DIR}/CommentLine.cpp
	${IO_SRC_DIR}/Prj.cpp
	${IO_SRC_DIR}/FileUtils.cpp
	${BASIC_SRC_DIR}/Word.cpp
)

add_executable(task_test
  ${DATABASE_TEST_SOURCES}
)

target_include_directories(task_test
  PRIVATE
    ${TASK_INCLUDE_DIR}
	${BASIC_INCLUDE_DIR}
	${IO_INCLUDE_DIR}
)

target_link_libraries(task_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(task_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
