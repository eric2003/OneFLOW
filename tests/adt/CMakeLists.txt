set(CMAKE_FOLDER ${CMAKE_FOLDER}/adt)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(ADT_INCLUDE_DIR ${OneFLOW_Root}/adt/include)
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(BASIC_SRC_DIR ${OneFLOW_Root}/basic/src)
set(GEOMETRY_INCLUDE_DIR ${OneFLOW_Root}/geometry/include)
set(GEOMETRY_SRC_DIR ${OneFLOW_Root}/geometry/src)
set(ADT_SRC_DIR ${OneFLOW_Root}/adt/src)

message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}" )
message(STATUS "ADT_INCLUDE_DIR=${ADT_INCLUDE_DIR}" )
message(STATUS "ADT_SRC_DIR=${ADT_SRC_DIR}" )

set(ADT_TEST_SOURCES
    adt_tree_test.cpp
    point_locator_test.cpp
	${GEOMETRY_SRC_DIR}/NodeMesh.cpp
	${GEOMETRY_SRC_DIR}/PointLocator.cpp
)

add_executable(adt_test
  ${ADT_TEST_SOURCES}
)

target_include_directories(adt_test
  PRIVATE
    ${ADT_INCLUDE_DIR}
	${ADT_SRC_DIR}
	${BASIC_INCLUDE_DIR}
	${GEOMETRY_INCLUDE_DIR}
)

target_link_libraries(adt_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(adt_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
