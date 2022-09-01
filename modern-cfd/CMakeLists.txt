cmake_minimum_required ( VERSION 3.20 )

project ( OneFLOW-1D )

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build directory: ${PROJECT_BINARY_DIR}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

list ( INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake" )

include(GNUInstallDirs)
include(util)

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

option ( PRJ_CUDA_SWITCH "Build Modern CFD(OneFLOW) with CUDA support" ON )
option ( PRJ_OPENMP_SWITCH "Build Modern CFD(OneFLOW) with OPENMP support" ON )
option ( PRJ_MPI_SWITCH "Build Modern CFD(OneFLOW) with MPI support" ON )

message ( STATUS  "PRJ_CUDA_SWITCH=${PRJ_CUDA_SWITCH}" )
message ( STATUS  "PRJ_OPENMP_SWITCH=${PRJ_OPENMP_SWITCH}" )
message ( STATUS  "PRJ_MPI_SWITCH=${PRJ_MPI_SWITCH}" )

if ( PRJ_MPI_SWITCH )
    find_package ( MPI )
    
    message ( STATUS  "MPI_FOUND=${MPI_FOUND}" )
    message ( STATUS  "MPI_CXX_INCLUDE_DIRS=${MPI_CXX_INCLUDE_DIRS}" )
    message ( STATUS  "MPI_LIBRARIES=${MPI_LIBRARIES}" )
    if ( MPI_FOUND )
        list ( APPEND PRJ_INCLUDE_DIRS ${MPI_CXX_INCLUDE_DIRS} )
        list ( APPEND PRJ_LIBRARIES ${MPI_LIBRARIES} )
        list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_MPI )
    endif ()
endif ()

if ( PRJ_CUDA_SWITCH )
    find_package ( CUDAToolkit )
    if ( CUDAToolkit_FOUND )
        enable_language( CUDA )
        set ( CMAKE_CUDA_STANDARD 20 )
        set ( CMAKE_CUDA_ARCHITECTURES "60;72;75" )
        list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_CUDA )
    endif ()
    message ( STATUS "CMAKE_CUDA_STANDARD=${CMAKE_CUDA_STANDARD}" )
endif ()

if ( PRJ_OPENMP_SWITCH )
    list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_OPENMP )
endif ()

if ( PRJ_OPENMP_SWITCH )
    find_package ( OpenMP )
    message ( STATUS  "OpenMP_FOUND=${OpenMP_FOUND}" )
    if ( OpenMP_FOUND )
        list ( APPEND PRJ_LIBRARIES OpenMP::OpenMP_CXX )
        list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_OPENMP )
        set ( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=${OpenMP_CXX_FLAGS}" )
    endif ()
endif ()

find_package ( Python3 COMPONENTS Interpreter Development NumPy )

message ( STATUS "Python3_FOUND = ${Python3_FOUND} " )
message ( STATUS "Python3_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS} " )
message ( STATUS "Python3_LIBRARIES = ${Python3_LIBRARIES} " )

#message ( STATUS "Python3_NumPy_FOUND = ${Python3_NumPy_FOUND} " )
#message ( STATUS "Python3_NumPy_INCLUDE_DIR = ${Python3_NumPy_INCLUDE_DIR} " )
#message ( STATUS "Python3_NumPy_INCLUDE_DIRS = ${Python3_NumPy_INCLUDE_DIRS} " )

list ( APPEND PRJ_INCLUDE_DIRS ${Python3_INCLUDE_DIRS} )
list ( APPEND PRJ_INCLUDE_DIRS ${Python3_NumPy_INCLUDE_DIRS} )
list ( APPEND PRJ_LIBRARIES ${Python3_LIBRARIES} )

get_directory_property( my_system_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} BUILDSYSTEM_TARGETS )
get_directory_property( my_import_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} IMPORTED_TARGETS )

message( STATUS "my_system_targets=${my_system_targets}" )
message( STATUS "my_import_targets=${my_import_targets}" )

message ( STATUS "MSVC=${MSVC}" )
if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

add_subdirectory ( src ${PROJECT_NAME} )
