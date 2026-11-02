# 根目录CMakeLists.txt
#cmake_minimum_required(VERSION 3.12)
cmake_minimum_required(VERSION 4.2.1)
project(FortranRegistry VERSION 1.0.0 LANGUAGES Fortran)

message(STATUS "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
message(STATUS "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")

# 设置Fortran标准
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 模块输出目录
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# 编译器选项
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
endif()

# 添加子目录
add_subdirectory(src)
add_subdirectory(tests)