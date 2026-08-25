# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)

project(OneFLOW_CFD LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  # 共享库(.so/.dylib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(SOURCES
	example/run_eno_weno.f90
)

add_subdirectory(src)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME}
    PRIVATE 
        cfd_lib
)

target_include_directories(${PROJECT_NAME} 
	PRIVATE 
		${CMAKE_Fortran_MODULE_DIRECTORY}
)
