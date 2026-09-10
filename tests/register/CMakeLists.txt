set(CMAKE_FOLDER ${CMAKE_FOLDER}/register)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
set(REGISTER_INCLUDE_DIR ${OneFLOW_Root}/register/include)
set(REGISTER_SRC_DIR ${OneFLOW_Root}/register/src)
set(IO_INCLUDE_DIR ${OneFLOW_Root}/io/include)
set(IO_SRC_DIR ${OneFLOW_Root}/io/src)

set(TEST_SOURCES
    MessageMapTest.cpp
	${REGISTER_SRC_DIR}/Message.cpp
	${IO_SRC_DIR}/TextFileParser.cpp
	${IO_SRC_DIR}/CommentLine.cpp
	${IO_SRC_DIR}/Prj.cpp
	${IO_SRC_DIR}/FileUtils.cpp
	${BASIC_SRC_DIR}/Word.cpp
	
)

add_executable(register_test
  ${TEST_SOURCES}
)

target_include_directories(register_test
  PRIVATE
    ${REGISTER_INCLUDE_DIR}
	${BASIC_INCLUDE_DIR}
	${IO_INCLUDE_DIR}
	
)

target_link_libraries(register_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(register_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
