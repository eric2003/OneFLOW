# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)

project(OneFLOW_CFD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 启用数学库
find_package(Threads REQUIRED)

# 源文件
set(SOURCES
    src/config.cpp
    src/mesh.cpp
    src/domain.cpp
    src/solution.cpp
    src/flux.cpp
    src/reconstructor.cpp
    src/cfd_solver.cpp
    src/main.cpp
)


add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(MSVC)
	target_compile_definitions( ${PROJECT_NAME} PRIVATE _USE_MATH_DEFINES )
    target_compile_options(${PROJECT_NAME} PRIVATE 
        "$<$<C_COMPILER_ID:MSVC>:/utf-8>"
        "$<$<CXX_COMPILER_ID:MSVC>:/utf-8>"
    )
	set_property(TARGET ${PROJECT_NAME} 
				 PROPERTY VS_STARTUP_PROJECT 
				 ${PROJECT_NAME}
	)
endif()

