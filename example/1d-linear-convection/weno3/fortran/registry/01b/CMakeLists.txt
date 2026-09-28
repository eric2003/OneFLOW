cmake_minimum_required(VERSION 4.2.1)
project(FortranCFD VERSION 1.0.0 LANGUAGES Fortran)

message(STATUS "Project: ${PROJECT_NAME} Version: ${PROJECT_VERSION}")
message(STATUS "Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Compiler Version: ${CMAKE_Fortran_COMPILER_VERSION}")

# Set Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Module output directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Compiler flags for Intel Fortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM|Intel")
    # Intel oneAPI/ifx compiler
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full /Od /traceback /check:all /warn:all /fpe:0")
    set(CMAKE_Fortran_FLAGS_RELEASE "/O3")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "/O2 /debug:full")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    # GNU gfortran compiler
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
endif()

# Set default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
    message(STATUS "Setting build type to: ${CMAKE_BUILD_TYPE}")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Create module directory
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

# Install directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install CACHE PATH "Installation prefix")

message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Module directory: ${CMAKE_Fortran_MODULE_DIRECTORY}")