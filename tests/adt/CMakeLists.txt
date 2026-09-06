set(CMAKE_FOLDER ${CMAKE_FOLDER}/adt)

set(OneFLOW_Root ${CMAKE_SOURCE_DIR}/codes )
set(ADT_INCLUDE_DIR ${OneFLOW_Root}/adt/include)
set(BASIC_INCLUDE_DIR ${OneFLOW_Root}/basic/include)
set(ADT_SRC_DIR ${OneFLOW_Root}/adt/src)

message(STATUS "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}" )
message(STATUS "ADT_INCLUDE_DIR=${ADT_INCLUDE_DIR}" )
message(STATUS "ADT_SRC_DIR=${ADT_SRC_DIR}" )

add_executable(adt_tree_test
  adt_tree_test.cpp
)

target_include_directories(adt_tree_test
  PRIVATE
    ${ADT_INCLUDE_DIR}
	${ADT_SRC_DIR}
	${BASIC_INCLUDE_DIR}
)

target_link_libraries(adt_tree_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(adt_tree_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
