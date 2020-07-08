project( Third_Party )

message( STATUS "CMAKE_PROJECT_NAME    = ${CMAKE_PROJECT_NAME}" )
message( STATUS "PROJECT_NAME = ${PROJECT_NAME}" )

GetAllSubDir( dirList ${CMAKE_CURRENT_SOURCE_DIR} )
message( STATUS "The dirList is ${dirList}" )

set(JSONCPP_DIR "JSONCpp")
set(TEST_JSON_DIR "TestJson")
add_subdirectory( jsoncpp ${JSONCPP_DIR} )
#add_subdirectory( testjson ${TEST_JSON_DIR} )
 