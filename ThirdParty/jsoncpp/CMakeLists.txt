project ( JSonCpp )

message ( STATUS "") 
message ( STATUS "CMAKE_CURRENT_LIST_FILE  = ${CMAKE_CURRENT_LIST_FILE}" )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}" )
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}" )

set ( JSonCpp_SRC_LIST jsoncpp.cpp )
set ( JSonCpp_HEAD_LIST
      json/json.h
      json/json-forwards.h
)

set ( BUILD_SHARED_JSON TRUE )
AppendGlobalValue ( THIRD_PROPERTY::INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} )

if ( BUILD_SHARED_JSON )
	add_library ( jsoncpp
		SHARED
		${JSonCpp_SRC_LIST}
		${JSonCpp_HEAD_LIST}
	)
	
	target_compile_definitions ( jsoncpp
		PRIVATE
			JSON_DLL_BUILD
	)
	
	AppendGlobalValue ( THIRD_PROPERTY::LIBRARIES   jsoncpp )
	AppendGlobalValue ( THIRD_PROPERTY::SHARED_LIBS jsoncpp )
	AppendGlobalValue ( THIRD_PROPERTY::DEFINE_SYMBOL JSON_DLL )
else ()
	add_library ( jsoncpp_static
		STATIC  
		${JSonCpp_SRC_LIST}
		${JSonCpp_HEAD_LIST}
	)
	AppendGlobalValue ( THIRD_PROPERTY::STATIC_LIBS jsoncpp_static )
	AppendGlobalValue ( THIRD_PROPERTY::LIBRARIES   jsoncpp_static )
endif ()

install ( DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/json" DESTINATION "include" )

 

