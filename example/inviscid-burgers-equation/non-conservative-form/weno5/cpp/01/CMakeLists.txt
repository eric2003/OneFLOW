cmake_minimum_required(VERSION 3.31)

project( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_23 )

message ( STATUS "EIGEN3_INCLUDE_DIR    = ${EIGEN3_INCLUDE_DIR}" )

message( STATUS "PRJ_INCLUDE_DIRS = ${PRJ_INCLUDE_DIRS}")

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

add_executable( ${PROJECT_NAME}
)

set( PROJECT_SOURCES
    main.cpp
)

message( STATUS "PROJECT_NAME = ${PROJECT_NAME}")
message( STATUS "PROJECT_SOURCES = ${PROJECT_SOURCES}")

target_sources( ${PROJECT_NAME}
  PRIVATE
	${PROJECT_SOURCES}
)

target_include_directories ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_INCLUDE_DIRS}
)

target_link_libraries( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_COMPILE_DEFINITIONS} 
)
