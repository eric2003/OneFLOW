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
option ( PRJ_FORTRAN_SWITCH "Build Modern CFD(OneFLOW) with CGNS support" OFF )
option ( PRJ_CGNS_SWITCH "Build Modern CFD(OneFLOW) with CGNS support" ON )
option ( PRJ_JSON_SWITCH "Build Modern CFD(OneFLOW) with JSON support" ON )

message ( STATUS  "PRJ_CUDA_SWITCH=${PRJ_CUDA_SWITCH}" )
message ( STATUS  "PRJ_OPENMP_SWITCH=${PRJ_OPENMP_SWITCH}" )
message ( STATUS  "PRJ_MPI_SWITCH=${PRJ_MPI_SWITCH}" )
message ( STATUS  "PRJ_FORTRAN_SWITCH=${PRJ_FORTRAN_SWITCH}" )
message ( STATUS  "PRJ_CGNS_SWITCH=${PRJ_CGNS_SWITCH}" )

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

if ( PRJ_CGNS_SWITCH )
    set ( CGNS_INCLUDE_DIRS $ENV{CGNS_INC} )
    set ( CGNS_LIBRARIES $ENV{CGNS_LIB_SHARED_NAME} )
    
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        set ( CGNS_INCLUDE_DIRS "C:/dev/cgns/include" CACHE PATH "path to CGNS headers" )
        set ( CGNS_LIBRARIES "C:/dev/cgns/lib/cgnsdll.lib" CACHE PATH "path to CGNS library" )
    endif()
    
    message(STATUS "CGNS_INCLUDE_DIRS = ${CGNS_INCLUDE_DIRS}")
    
    list ( APPEND PRJ_LIBRARIES ${CGNS_LIBRARIES} )
    list ( APPEND PRJ_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} )
    list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_CGNS )
endif ()

if ( PRJ_JSON_SWITCH )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        set ( debug_dir $<$<CONFIG:Debug>:/debug> )
        set ( JSON_INCLUDE_DIRS "C:/dev${debug_dir}/jsoncpp/include" CACHE PATH "path to JSON headers" )
        set ( JSON_LIBRARIES "C:/dev${debug_dir}/jsoncpp/lib/jsoncpp.lib" CACHE PATH "path to JSON library" )
    endif()
    
    list ( APPEND PRJ_LIBRARIES ${JSON_LIBRARIES} )
    list ( APPEND PRJ_INCLUDE_DIRS ${JSON_INCLUDE_DIRS} )
    list ( APPEND PRJ_COMPILE_DEFINITIONS PRJ_ENABLE_JSON )
endif ()

get_directory_property( my_system_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} BUILDSYSTEM_TARGETS )
get_directory_property( my_import_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} IMPORTED_TARGETS )

message( STATUS "my_system_targets=${my_system_targets}" )
message( STATUS "my_import_targets=${my_import_targets}" )

message ( STATUS "MSVC=${MSVC}" )
if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

add_subdirectory ( src ${PROJECT_NAME} )
