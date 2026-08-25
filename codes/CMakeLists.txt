message ( STATUS "") 
message ( STATUS "CMAKE_CURRENT_LIST_FILE = ${CMAKE_CURRENT_LIST_FILE}." )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}." )

set ( PRJ_INCLUDE_DIRS )
set ( PRJ_LIBRARIES )
set ( PRJ_LINK_DIRS )
set ( PRJ_COMPILE_DEF )
set ( PRJ_COMPILE_FEATURES )
set ( USER_LIBS )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_17 )

list ( APPEND PRJ_COMPILE_DEF _HAS_STD_BYTE=0 )

set ( oneAPI_LIB_DIR )
set ( LIB_COMPILER_TYPE VS2022 )
#set ( LIB_COMPILER_TYPE oneAPI2022 )
set ( LIB_TYPE STATIC )
if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
    #Cancel the warning of fopen in C language
    list ( APPEND PRJ_COMPILE_DEF _CRT_SECURE_NO_WARNINGS )
    set ( oneAPI_LIB_DIR "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/compiler/lib/intel64/" )
    list ( APPEND PRJ_LINK_DIRS ${oneAPI_LIB_DIR} )
    #list ( APPEND PRJ_COMPILE_DEF USE_METIS_DLL )
endif()

set ( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code" )
set ( AUTO_ACC TRUE CACHE BOOL "auto acc or not")

set ( ONEFLOW_ACCEL_BACKEND "AUTO" CACHE STRING
      "Default OneFLOW compute backend (AUTO, CPU, HIP, CUDA, or KOKKOS)" )
set_property ( CACHE ONEFLOW_ACCEL_BACKEND PROPERTY STRINGS AUTO CPU HIP CUDA KOKKOS )
option ( ONEFLOW_ENABLE_MULTI_DEVICE
         "Keep rank/device mapping information in the backend runtime" ON )
option ( ONEFLOW_ENABLE_HIP
         "Explicitly build the native HIP backend" OFF )
option ( ONEFLOW_AUTO_ENABLE_HIP
         "Build HIP automatically when compiler, package, and device architecture are detected" ON )
option ( ONEFLOW_ENABLE_CUDA
         "Build the legacy CUDA sources explicitly" OFF )

string ( TOUPPER "${ONEFLOW_ACCEL_BACKEND}" ONEFLOW_ACCEL_BACKEND_UPPER )
set ( ONEFLOW_ACCEL_BACKEND_VALUES AUTO CPU HIP CUDA KOKKOS )
if ( NOT ONEFLOW_ACCEL_BACKEND_UPPER IN_LIST ONEFLOW_ACCEL_BACKEND_VALUES )
    message ( FATAL_ERROR
        "Unsupported ONEFLOW_ACCEL_BACKEND=${ONEFLOW_ACCEL_BACKEND}. "
        "Expected one of: AUTO, CPU, HIP, CUDA, KOKKOS." )
endif()
if ( ONEFLOW_ACCEL_BACKEND_UPPER STREQUAL "CUDA"
     AND NOT ONEFLOW_ENABLE_CUDA )
    message ( FATAL_ERROR
        "CUDA is a reserved extension point. Set ONEFLOW_ENABLE_CUDA=ON "
        "only when the legacy CUDA source path is intentionally enabled." )
endif()
if ( ONEFLOW_ACCEL_BACKEND_UPPER STREQUAL "KOKKOS" )
    message ( FATAL_ERROR
        "KOKKOS is reserved but is not enabled in the current OneFLOW port." )
endif()

list ( APPEND PRJ_COMPILE_DEF
       ONEFLOW_DEFAULT_ACCEL_BACKEND="${ONEFLOW_ACCEL_BACKEND_UPPER}" )
if ( ONEFLOW_ENABLE_MULTI_DEVICE )
    list ( APPEND PRJ_COMPILE_DEF ONEFLOW_ENABLE_MULTI_DEVICE )
endif()

set ( ONEFLOW_BUILD_HIP OFF )
set ( ONEFLOW_HIP_EXPLICIT OFF )
if ( ONEFLOW_ENABLE_HIP OR ONEFLOW_ACCEL_BACKEND_UPPER STREQUAL "HIP" )
    set ( ONEFLOW_HIP_EXPLICIT ON )
endif()

if ( ONEFLOW_HIP_EXPLICIT
     OR ( ONEFLOW_AUTO_ENABLE_HIP
          AND ONEFLOW_ACCEL_BACKEND_UPPER STREQUAL "AUTO" ) )
    DetectHIPArchitectures ( ONEFLOW_DETECTED_HIP_ARCHITECTURES )
    if ( ONEFLOW_DETECTED_HIP_ARCHITECTURES )
        set ( CMAKE_HIP_ARCHITECTURES "${ONEFLOW_DETECTED_HIP_ARCHITECTURES}"
              CACHE STRING "HIP/DCU architecture list" FORCE )
    endif()

    find_package ( hip CONFIG QUIET )
    include ( CheckLanguage )
    check_language ( HIP )
    if ( hip_FOUND AND CMAKE_HIP_COMPILER
         AND CMAKE_HIP_ARCHITECTURES )
        set ( ONEFLOW_BUILD_HIP ON )
    elseif ( ONEFLOW_HIP_EXPLICIT )
        message ( FATAL_ERROR
            "HIP was explicitly requested, but the HIP package, compiler, or "
            "device architecture was not detected. Set CMAKE_HIP_COMPILER and "
            "CMAKE_HIP_ARCHITECTURES explicitly when cross-compiling." )
    endif()
endif()

message ( STATUS "ONEFLOW_ACCEL_BACKEND=${ONEFLOW_ACCEL_BACKEND_UPPER}" )
message ( STATUS "ONEFLOW_ENABLE_MULTI_DEVICE=${ONEFLOW_ENABLE_MULTI_DEVICE}" )
message ( STATUS "ONEFLOW_ENABLE_HIP=${ONEFLOW_ENABLE_HIP}" )
message ( STATUS "ONEFLOW_AUTO_ENABLE_HIP=${ONEFLOW_AUTO_ENABLE_HIP}" )
message ( STATUS "ONEFLOW_BUILD_HIP=${ONEFLOW_BUILD_HIP}" )
message ( STATUS "ONEFLOW_ENABLE_CUDA=${ONEFLOW_ENABLE_CUDA}" )
if ( ONEFLOW_BUILD_HIP )
    message ( STATUS "CMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES}" )
    enable_language ( HIP )
    list ( APPEND PRJ_COMPILE_DEF ONEFLOW_ENABLE_HIP )
    list ( APPEND PRJ_INCLUDE_DIRS ${hip_INCLUDE_DIRS} )
    list ( APPEND PRJ_LIBRARIES hip::host )

    GetHIPFiles(
        ${CMAKE_CURRENT_SOURCE_DIR}
        HIP_HEADER_LIST
        HIP_SRC_LIST
        HIP_LOCAL_INCLUDE_DIRS )

    add_subdirectory ( accel )
endif()

if ( ONEFLOW_ENABLE_CUDA )
    find_package ( CUDAToolkit REQUIRED )
else()
    set ( CUDAToolkit_FOUND OFF )
endif()
find_package ( OpenACC )

message ( STATUS "OpenACC_FOUND      = ${OpenACC_FOUND}" )
message ( STATUS "OpenACC_C_FOUND    = ${OpenACC_C_FOUND}" )
message ( STATUS "OpenACC_CXX_FOUND  = ${OpenACC_CXX_FOUND}" )
message ( STATUS "OpenACC_C_FLAGS    = ${OpenACC_C_FLAGS}" )
message ( STATUS "OpenACC_CXX_FLAGS  = ${OpenACC_CXX_FLAGS}" )

#set ( GPUFLAGS "-gpu=managed -Minfo=all,ccff -Mneginfo" )
#set ( GPUFLAGS "-gpu=managed -Minfo=all,ccff" )
set ( GPUFLAGS "-gpu=managed -Minfo=accel" )

if ( OpenACC_CXX_FOUND AND AUTO_ACC )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_OPENACC )
    message ( STATUS "ENV{NVCC_COMPILER}  = $ENV{NVCC_COMPILER}" )
    if ( $ENV{NVCC_COMPILER} MATCHES "NVCC_COMPILER" )
        #set( ACC_FLAG "-fast -Minfo=accel -gpu=cc60,cc70,cc75,cc80" )
        set( ACC_FLAG "-fast ${GPUFLAGS}" )
        list ( APPEND PRJ_COMPILE_DEF NVCC_COMPILER_ENABLED )
    else()
        #set( ACC_FLAG "-Ofast" )
        #set( ACC_FLAG "-Ofast ${GPUFLAGS}" )
        #set( ACC_FLAG "-Ofast -foffload=disable" )
        set( ACC_FLAG "-Ofast -foffload=disable" )
    endif()
    set ( OpenACC_CXX_FLAGS "${OpenACC_CXX_FLAGS} ${ACC_FLAG}" )
endif ()

message ( STATUS "OpenACC_CXX_FLAGS now = ${OpenACC_CXX_FLAGS}" )

find_package ( OpenMP )
if ( OpenMP_FOUND )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_OPENMP )
    message ( STATUS "OpenMP_CXX_INCLUDE_DIRS  = ${OpenMP_CXX_INCLUDE_DIRS}" )
    message ( STATUS "OpenMP_CXX_LIBRARIES  = ${OpenMP_CXX_LIBRARIES}" )
    message ( STATUS "OpenMP_CXX_FLAGS  = ${OpenMP_CXX_FLAGS}" )
    list ( APPEND PRJ_LIBRARIES ${OpenMP_CXX_LIBRARIES} )
endif ()

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${OpenACC_CXX_FLAGS}" )

message ( STATUS "CMAKE_CXX_FLAGS now = ${CMAKE_CXX_FLAGS}" )

if ( CUDAToolkit_FOUND )
    enable_language( CUDA )
    set ( CMAKE_CUDA_STANDARD 14 )
    set ( CMAKE_CUDA_ARCHITECTURES "60;72;75" )
    set ( CUDA_INCLUDE_DIRS )
    set ( CUDA_PRJ_COMPILE_DEF )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_CUDA )
    list ( APPEND CUDA_PRJ_COMPILE_DEF ${PRJ_COMPILE_DEF} )
    
    GetCUDAFiles( ${CMAKE_CURRENT_SOURCE_DIR} CUDA_HEADER_LIST CUDA_SRC_LIST CUDA_LOCAL_INCLUDE_DIRS )
    
    list ( APPEND CUDA_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} )
    list ( APPEND CUDA_INCLUDE_DIRS ${CUDA_LOCAL_INCLUDE_DIRS} )
    
    set ( CUDA_LIB "CUDA_LIB" )
    add_library( ${CUDA_LIB} STATIC 
        ${CUDA_SRC_LIST}
    )
    target_include_directories ( ${CUDA_LIB}
        PRIVATE 
            ${CUDA_INCLUDE_DIRS}
    )
    
    target_compile_definitions ( ${CUDA_LIB}
        PRIVATE
           ${CUDA_PRJ_COMPILE_DEF} 
    )
    
    list ( APPEND PRJ_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS} )
    list ( APPEND PRJ_LIBRARIES ${CUDA_LIB} )
    list ( APPEND USER_LIBS ${CUDA_LIB} )
else()
    message ( STATUS "No CUDA support" )
endif()

message ( STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}" )
message ( STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}" )
message ( STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}" )
message ( STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}" )
message ( STATUS "The CMAKE_EXE_LINKER_FLAGS is ${CMAKE_EXE_LINKER_FLAGS}" )

message ( STATUS "The CFLAGS is ${CFLAGS}" )
message ( STATUS "The CXXFLAGS is ${CXXFLAGS}" )
message ( STATUS "The LDFLAGS is ${LDFLAGS}" )

if ( BIT64_FLAG )
    set ( ARCH "x64" )
else()
    set ( ARCH "x86" )
endif()

#MPI library settings
option ( MPI_ENABLE "ON for using MPI library" ON )

if ( MPI_ENABLE )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_MPI )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        set ( MPI_INCLUDE_DIRS "C:/Program Files (x86)/Microsoft SDKs/MPI/include" CACHE PATH "path to mpi headers" )
        set ( MPI_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${ARCH}/msmpi.lib" CACHE PATH "path to mpi library" )
    else()
        set ( MPI_INCLUDE_DIRS "$ENV{MPI_HOME_INC}" CACHE PATH "path to mpi headers" )
        set ( MPI_LIBRARIES "$ENV{MPI_HOME_LIB}" CACHE PATH "path to mpi library" )
    endif()
    
    message ( STATUS "The MPI_LIBRARIES is ${MPI_LIBRARIES}" )
    
    list ( APPEND PRJ_LIBRARIES ${MPI_LIBRARIES} )
    list ( APPEND PRJ_INCLUDE_DIRS ${MPI_INCLUDE_DIRS} )
endif()


#METIS library settings
option ( METIS_ENABLE "ON for using METIS partitioning library" ON )

if ( METIS_ENABLE )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_METIS )
    
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set ( METIS_ROOT "C:/dev/METIS/METIS-${LIB_COMPILER_TYPE}-STATIC" )
        set ( METIS_INCLUDE_DIRS "${METIS_ROOT}/include"  CACHE PATH "path to METIS headers" )
        set ( METIS_LIBRARIES "${METIS_ROOT}/lib/metis.lib" CACHE PATH "path to METIS library" )
    else()
        set ( METIS_INCLUDE_DIRS $ENV{METIS_HOME_INC} )
        set ( METIS_LIBRARIES $ENV{METIS_HOME_LIB} )
    endif()

    list ( APPEND PRJ_LIBRARIES ${METIS_LIBRARIES} )
    list ( APPEND PRJ_INCLUDE_DIRS ${METIS_INCLUDE_DIRS} )
endif()

#CGNS library settings
option ( CGNS_ENABLE "ON for using CGNS library" ON )

if ( CGNS_ENABLE )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_CGNS )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        set ( CGNS_ROOT "C:/dev/cgns/4.4.0" )
        set ( CGNS_LIBNAME "cgnsdll.lib" )
        set ( CGNS_INCLUDE_DIRS "${CGNS_ROOT}/include" CACHE PATH "path to CGNS headers" )
        set ( CGNS_LIBRARIES "${CGNS_ROOT}/lib/${CGNS_LIBNAME}" CACHE PATH "path to CGNS library" )	
    else()
        set ( CGNS_INCLUDE_DIRS $ENV{CGNS_HOME_INC} CACHE PATH "path to CGNS headers" )
        set ( CGNS_LIBRARIES $ENV{CGNS_HOME_LIB} CACHE PATH "path to CGNS library" )
    endif()
    message(STATUS "The CGNS_INCLUDE_DIRS is ${CGNS_INCLUDE_DIRS}")
    
    list ( APPEND PRJ_LIBRARIES ${CGNS_LIBRARIES} )
    list ( APPEND PRJ_INCLUDE_DIRS ${CGNS_INCLUDE_DIRS} )
endif()

set ( JsonPrjDir "${PRJ_ROOT_DIR}/${THIRD_PARTY_DIR}/jsoncpp" )
set ( JsonHeaderDir ${JsonPrjDir} )

ConstructSolutionDirTree( ${CMAKE_CURRENT_SOURCE_DIR} PRJ_HEADER_LIST PRJ_SRC_LIST PRJ_LOCAL_INCLUDE_DIRS )
list ( FILTER PRJ_SRC_LIST EXCLUDE REGEX "/codes/accel/src/HipSmoke\\.cpp$" )
if ( ONEFLOW_BUILD_HIP )
    list ( APPEND PRJ_HEADER_LIST ${HIP_HEADER_LIST} )
    list ( APPEND PRJ_SRC_LIST ${HIP_SRC_LIST} )
    list ( APPEND PRJ_INCLUDE_DIRS ${HIP_LOCAL_INCLUDE_DIRS} )
endif()
list ( APPEND PRJ_INCLUDE_DIRS ${PRJ_LOCAL_INCLUDE_DIRS} )

if ( NOT CGNS_ENABLE )
    list ( FILTER PRJ_SRC_LIST EXCLUDE REGEX "/codes/cgns/" )
    list ( FILTER PRJ_HEADER_LIST EXCLUDE REGEX "/codes/cgns/" )
    list ( FILTER PRJ_INCLUDE_DIRS EXCLUDE REGEX "/codes/cgns/" )
endif()

GetGlobalValue ( THIRD_PROPERTY::INCLUDE_DIRS THIRD_INCLUDE_DIRS )

list ( APPEND PRJ_INCLUDE_DIRS ${THIRD_INCLUDE_DIRS} )

message( STATUS "THIRD_INCLUDE_DIRS = ${THIRD_INCLUDE_DIRS}" )

GetGlobalValue ( THIRD_PROPERTY::LIBRARIES THIRD_LIBRARIES )
GetGlobalValue ( THIRD_PROPERTY::DEFINE_SYMBOL THIRD_DEFINE_SYMBOL )

list ( APPEND PRJ_LIBRARIES ${THIRD_LIBRARIES} )
list ( APPEND PRJ_COMPILE_DEF ${THIRD_DEFINE_SYMBOL} )

if ( BIT64_FLAG )
    list ( APPEND PRJ_COMPILE_DEF X64 )
endif( BIT64_FLAG )

if ( WIN32 )
    list ( APPEND PRJ_COMPILE_DEF HX_PARALLEL )
    list ( APPEND PRJ_COMPILE_DEF HX_CGNS )
    list ( APPEND PRJ_COMPILE_DEF USE_DLL )
elseif ( MPI_ENABLE )
    list ( APPEND PRJ_COMPILE_DEF HX_PARALLEL )
    list ( APPEND PRJ_COMPILE_DEF HX_CGNS )
endif ()

add_executable ( ${PROJECT_NAME} ${PRJ_SRC_LIST} ${PRJ_HEADER_LIST} )

target_include_directories ( ${PROJECT_NAME}
    PRIVATE 
        ${PRJ_INCLUDE_DIRS}
)

target_link_directories ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LINK_DIRS}
)

target_link_libraries ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
    PRIVATE
       ${PRJ_COMPILE_DEF} 
)

set ( PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/ )
set ( PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/ )

install ( TARGETS ${PROJECT_NAME} ${USER_LIBS} ${THIRD_LIBRARIES}
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install ( DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system" )
install ( DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc" )

message ( STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---" )
