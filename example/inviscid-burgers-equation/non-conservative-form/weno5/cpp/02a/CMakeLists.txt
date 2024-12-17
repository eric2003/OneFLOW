cmake_minimum_required(VERSION 3.31)

project( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_23 )

message ( STATUS "EIGEN3_INCLUDE_DIR    = ${EIGEN3_INCLUDE_DIR}" )

set ( CGNS_INCLUDE_DIRS $ENV{CGNS_INC} )
set ( CGNS_LIBRARIES $ENV{CGNS_LIB_SHARED_NAME} )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    set ( CGNS_ROOT "C:/dev/cgns/4.4.0" )
    set ( CGNS_LIBNAME "cgnsdll.lib" )

    set ( CGNS_INCLUDE_DIRS "${CGNS_ROOT}/include" CACHE PATH "path to CGNS headers" )
    set ( CGNS_LIBRARIES "${CGNS_ROOT}/lib/${CGNS_LIBNAME}" CACHE PATH "path to CGNS library" )
endif()

list ( APPEND PRJ_LIBRARIES ${CGNS_LIBRARIES} )
list ( APPEND PRJ_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} )

if ( WIN32 )
    list ( APPEND PRJ_COMPILE_DEFINITIONS USE_DLL )
endif ()

#list ( APPEND PRJ_COMPILE_DEFINITIONS HX_PARALLEL )

message( STATUS "PRJ_INCLUDE_DIRS = ${PRJ_INCLUDE_DIRS}")

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

add_executable( ${PROJECT_NAME}
)


set( PROJECT_SOURCES
    main.cpp
	Field.h Field.cpp
	Global.h Global.cpp
	Grid.h Grid.cpp
	LogFile.h LogFile.cpp
	Parallel.h Parallel.cpp
	Post.h Post.cpp
	Solver.h Solver.cpp
	Weno.h Weno.cpp
    ZoneState.h ZoneState.cpp
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
