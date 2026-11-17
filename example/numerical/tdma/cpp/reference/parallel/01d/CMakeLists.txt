cmake_minimum_required(VERSION 3.31)
project(testprj)

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_23 )

find_package ( MPI )

message ( STATUS  "MPI_FOUND=${MPI_FOUND}" )
message ( STATUS  "MPI_CXX_INCLUDE_DIRS=${MPI_CXX_INCLUDE_DIRS}" )
message ( STATUS  "MPI_LIBRARIES=${MPI_LIBRARIES}" )
if ( MPI_FOUND )
	list ( APPEND PRJ_INCLUDE_DIRS ${MPI_CXX_INCLUDE_DIRS} )
	list ( APPEND PRJ_LIBRARIES ${MPI_LIBRARIES} )
endif ()

list ( APPEND PRJ_LIBRARIES MPI::MPI_C )

add_executable( ${PROJECT_NAME} )

target_sources( ${PROJECT_NAME}
    PRIVATE
        main.cpp
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