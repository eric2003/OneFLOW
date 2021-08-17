project ( Third_Party )

message ( STATUS "CMAKE_PROJECT_NAME    = ${CMAKE_PROJECT_NAME}" )
message ( STATUS "PROJECT_NAME = ${PROJECT_NAME}" )

set ( JSONCPP_DIR "JSONCpp" )
set ( TEST_JSON_DIR "TestJson" )
add_subdirectory ( jsoncpp ${JSONCPP_DIR} )
 