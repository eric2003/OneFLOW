set(CMAKE_FOLDER ${CMAKE_FOLDER}/database)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(DATABASE_INCLUDE_DIR ${OneFLOW_Root}/database/include)
set(DATABASE_SRC_DIR ${OneFLOW_Root}/database/src)
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
#set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
#set(GEOMETRY_INCLUDE_DIR ${OneFLOW_Root}/geometry/include)
#set(GEOMETRY_SRC_DIR ${OneFLOW_Root}/geometry/src)
#set(ADT_SRC_DIR ${OneFLOW_Root}/adt/src)
#
#message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}" )
#message(STATUS "ADT_INCLUDE_DIR=${ADT_INCLUDE_DIR}" )
#message(STATUS "ADT_SRC_DIR=${ADT_SRC_DIR}" )

set(DATABASE_TEST_SOURCES
    database_para_test.cpp
)

add_executable(database_test
  ${DATABASE_TEST_SOURCES}
)

target_include_directories(database_test
  PRIVATE
    ${DATABASE_INCLUDE_DIR}
	${BASIC_INCLUDE_DIR}
)

target_link_libraries(database_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(database_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
