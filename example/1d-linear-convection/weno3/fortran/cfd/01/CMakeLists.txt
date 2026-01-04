# CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)  # 更高版本更好支持Fortran

project(OneFLOW_CFD LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)  # 共享库(.so/.dylib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 按依赖顺序显式列出源文件（解决 Intel + VS 编译顺序问题）
set(SOURCES
    cfd_solver.f90
)

add_executable(oneflow_cfd ${SOURCES})

# 包含模块目录（解决 .mod 未找到）
target_include_directories(oneflow_cfd PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Intel 特定（可选，链接数学库）
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_link_libraries(oneflow_cfd m)
endif()