cmake_minimum_required(VERSION 3.20)

project ( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )

enable_language(Fortran)

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

add_executable( ${PROJECT_NAME}
    main.f90
)

target_compile_features ( ${PROJECT_NAME} 
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)