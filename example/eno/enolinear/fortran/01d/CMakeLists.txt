cmake_minimum_required(VERSION 3.31)

project ( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )

enable_language(Fortran)

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

add_executable( ${PROJECT_NAME}
    enolinear.f90
)

target_compile_features ( ${PROJECT_NAME} 
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)