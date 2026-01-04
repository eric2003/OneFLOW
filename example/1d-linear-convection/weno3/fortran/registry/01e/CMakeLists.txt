# 根目录CMakeLists.txt
cmake_minimum_required(VERSION 4.2.1)
project(FortranCFD VERSION 1.0.0 LANGUAGES Fortran)

message(STATUS "Project: ${PROJECT_NAME} Version: ${PROJECT_VERSION}")
message(STATUS "Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Compiler Version: ${CMAKE_Fortran_COMPILER_VERSION}")

# 设置Fortran标准
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# 模块输出目录
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# 编译器选项
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full /Od /traceback /check:all /warn:all /fpe:0")
    set(CMAKE_Fortran_FLAGS_RELEASE "/O3")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "/O2 /debug:full")
endif()

# 设置默认构建类型
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
    message(STATUS "Setting default build type to: ${CMAKE_BUILD_TYPE}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# 确保模块目录存在
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# 添加子目录
add_subdirectory(src)
add_subdirectory(tests)

# 安装目录
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Installation prefix")

message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")