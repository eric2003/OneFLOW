project ( JSonCpp )

message( STATUS "CMAKE_PROJECT_NAME = ${CMAKE_PROJECT_NAME}" )
message( STATUS "PROJECT_NAME       = ${PROJECT_NAME}" )

set( JSon_EXE "${CMAKE_PROJECT_NAME}" )
message( STATUS "The JSon_EXE is ${JSon_EXE}.")

set ( JSonCpp_SRC_LIST jsoncpp.cpp )
set ( JSonCpp_HEAD_LIST
      json/json.h
      json/json-forwards.h
)

#Generate the static library from the library sources
add_library ( ${PROJECT_NAME} STATIC  ${JSonCpp_SRC_LIST} ${JSonCpp_HEAD_LIST} )
 

