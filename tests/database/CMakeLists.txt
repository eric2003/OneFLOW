set(CMAKE_FOLDER ${CMAKE_FOLDER}/database)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(DATABASE_INCLUDE_DIR ${OneFLOW_Root}/database/include)
set(DATABASE_SRC_DIR ${OneFLOW_Root}/database/src)
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
set(ARRAY_INCLUDE_DIR ${OneFLOW_Root}/array/include)
set(PARALLEL_INCLUDE_DIR ${OneFLOW_Root}/parallel/include)
set(PARALLEL_SRC_DIR ${OneFLOW_Root}/parallel/src)
set(ZONE_INCLUDE_DIR ${OneFLOW_Root}/zone/include)
set(ZONE_SRC_DIR ${OneFLOW_Root}/zone/src)
set(GLOBAL_INCLUDE_DIR ${OneFLOW_Root}/global/include)

set(DATABASE_TEST_SOURCES
    database_para_test.cpp
	database_field_test.cpp
	${DATABASE_SRC_DIR}/DataBase.cpp
	${DATABASE_SRC_DIR}/DataPara.cpp
	${DATABASE_SRC_DIR}/DataBook.cpp
	${DATABASE_SRC_DIR}/DataField.cpp
	${DATABASE_SRC_DIR}/DataBaseType.cpp
	${DATABASE_SRC_DIR}/DataPage.cpp
	${DATABASE_SRC_DIR}/DataBaseIO.cpp
	${BASIC_SRC_DIR}/Stop.cpp
	${PARALLEL_SRC_DIR}/BasicParallel.cpp
	${PARALLEL_SRC_DIR}/Parallel.cpp
	${ZONE_SRC_DIR}/ZoneState.cpp
)

add_executable(database_test
  ${DATABASE_TEST_SOURCES}
)

target_include_directories(database_test
  PRIVATE
    ${DATABASE_INCLUDE_DIR}
	${BASIC_INCLUDE_DIR}
	${ARRAY_INCLUDE_DIR}
	${PARALLEL_INCLUDE_DIR}
	${ZONE_INCLUDE_DIR}
	${GLOBAL_INCLUDE_DIR}
)

target_link_libraries(database_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(database_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
