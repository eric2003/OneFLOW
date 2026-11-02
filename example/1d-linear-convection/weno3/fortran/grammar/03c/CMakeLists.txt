cmake_minimum_required(VERSION 3.20)

project ( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )

enable_language(Fortran)

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )
list ( APPEND PRJ_LIBRARIES mylib )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  # 共享库(.so/.dylib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(src)

add_executable( ${PROJECT_NAME}
    main.f90
)

target_include_directories( ${PROJECT_NAME} 
	PRIVATE 
		${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_compile_features ( ${PROJECT_NAME} 
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)

target_link_libraries( ${PROJECT_NAME}
	PRIVATE 
		${PRJ_LIBRARIES} 
)




