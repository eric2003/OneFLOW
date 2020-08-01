#. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code")

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if (CMAKE_GENERATOR MATCHES "Visual Studio 16 2019" )
        set(VS_COMPILER vs2019)
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 15 2017" )
        set(VS_COMPILER vs2017)
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 14 2015" )
        set(VS_COMPILER vs2015)
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 12 2013" )
        set(VS_COMPILER vs2013)
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 11 2012" )
        set(VS_COMPILER vs2012)
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 10 2010" )
        set(VS_COMPILER vs2010)
    else()
        set(VS_COMPILER vs2019)
    endif(CMAKE_GENERATOR MATCHES "Visual Studio 16 2019" )
    message(STATUS "The current VS_COMPILER is ${VS_COMPILER}")
endif()

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    #set(CMAKE_CXX_FLAGS "-std=c++11" CACHE STRING "" FORCE)
    option(USE_PREBUILT_LIB "ON for using prebuilt third party libraries" OFF)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "VS_COMPILER is ${VS_COMPILER}")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)  #取消c语言中fopen等函数的警告
    option(USE_PREBUILT_LIB "ON for using prebuilt third party libraries" ON) #只在windows下使用已编译好的第三方库
endif()

message(STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}")
message(STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message(STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message(STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

message(STATUS "The CFLAGS is ${CFLAGS}")
message(STATUS "The CXXFLAGS is ${CXXFLAGS}")
message(STATUS "The LDFLAGS is ${LDFLAGS}")

message(STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message(STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")
message(STATUS "The CMAKE_EXE_LINKER_FLAGS is ${CMAKE_EXE_LINKER_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -lrt")
Set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -ldl")
message(STATUS "The CMAKE_EXE_LINKER_FLAGS is ${CMAKE_EXE_LINKER_FLAGS}")

if( BIT64_FLAG )
    set(ARCH "x64")
else()
    set(ARCH "x86")
endif()

#使用已经编译好的库
if( USE_PREBUILT_LIB )
    set( ThirdPartyLib "ThirdPartyLib" )
    set(THIRD_LIB_BASE_DIR "${CMAKE_SOURCE_DIR}/${ThirdPartyLib}")
    set(THIRD_LIB_DIR "${CMAKE_SOURCE_DIR}/${ThirdPartyLib}/${VS_COMPILER}/${ARCH}")
endif()

#mpi库设置(不需要区别vs版本)
option(MPI_ENABLE "ON for using MPI library" ON)

if ( MPI_ENABLE )
    add_definitions(-DENABLE_MPI)
    if ( MPI_DIR )
    else()
        if ( WIN32 )
            if ( AUTO_CI_TEST )
                set( MPI_DIR "C:/Program Files (x86)/Microsoft SDKs/MPI" )
            else()
                set( MPI_DIR "${THIRD_LIB_BASE_DIR}/msmpi" )
            endif()
        else()
            set( MPI_DIR "/usr/local/mpich" )
            if ( AUTO_CI_TEST )
                set( MPI_DIR "/home/travis/build/eric2003/OneFLOW/openmpi" )
            endif()
        endif()
    endif()

    if  ( USE_PREBUILT_LIB AND WIN32 )
        set(MPI_HEADER_DIR "${MPI_DIR}/include" CACHE PATH "path to mpi headers")
        set(MPI_LIBRARY_DIR "${MPI_DIR}/lib/${ARCH}" CACHE PATH "path to mpi library")
    else()
        set(MPI_HEADER_DIR "${MPI_DIR}/include" CACHE PATH "path to mpi headers")
        set(MPI_LIBRARY_DIR "${MPI_DIR}/lib" CACHE PATH "path to mpi library")
    endif()

    message(STATUS "AUTO_CI_TEST value is ${AUTO_CI_TEST}")

    if( AUTO_CI_TEST )
        # determine if MPI_DIR has been set
        #-DMPI_DIR=$TRAVIS_BUILD_DIR/openmpi
        if (MPI_DIR)
            message (STATUS "Looking for MPI here: ${MPI_DIR}")
            set (CMAKE_PREFIX_PATH ${MPI_DIR})
        endif ()

        message(STATUS "MPI_HEADER_DIR is ${MPI_HEADER_DIR}")
        message(STATUS "MPI_LIBRARY_DIR is ${MPI_LIBRARY_DIR}")

        #include_directories (SYSTEM "${MPI_INCLUDE_PATH}")
        #include_directories ("${MPI_HEADER_DIR}")

        #message( STATUS "MPI_DIR is ${MPI_DIR}" )
        #message( STATUS "MPI_INCLUDE_PATH is ${MPI_INCLUDE_PATH}" )
        #message( STATUS "MPI_C_LIBRARIES is ${MPI_C_LIBRARIES}" )

        #message( STATUS "MPI_HEADER_DIR is ${MPI_HEADER_DIR}" )
        #message( STATUS "MPI_LIBRARY_DIR is ${MPI_LIBRARY_DIR}" )
    endif ()

    if( WIN32 )
        set(MPI_LIBRARY_NAME "msmpi" CACHE STRING "name of mpi library")
    else()
        set(MPI_LIBRARY_NAME "mpi" CACHE STRING "name of mpi library")
        if( AUTO_CI_TEST )
            set(MPI_LIBRARY_NAME "libmpi.so" CACHE STRING "name of mpi library" FORCE)
        endif()
    endif()

    include_directories(${MPI_HEADER_DIR})  #头文件目录
    link_directories(${MPI_LIBRARY_DIR})    #库目录
    link_libraries(${MPI_LIBRARY_NAME})     #库名称

    message(STATUS "MPI_HEADER_DIR is ${MPI_HEADER_DIR}")
    message(STATUS "MPI_LIBRARY_DIR is ${MPI_LIBRARY_DIR}")
    message(STATUS "MPI_LIBRARY_NAME is ${MPI_LIBRARY_NAME}")
endif()

#metis库设置
option(METIS_ENABLE "ON for using METIS partitioning library" ON)

if(METIS_ENABLE)
    add_definitions(-DENABLE_METIS)
    
    if(USE_PREBUILT_LIB)
        set(METIS_DIR "${THIRD_LIB_DIR}/metis")
        if( AUTO_CI_TEST )
            set(METIS_DIR "C:/Program Files/metis")
        endif()
    else()
        set(METIS_DIR "/usr/local/metis")
    endif()

    set(METIS_HEADER_DIR "${METIS_DIR}/include" CACHE PATH "path to metis headers")
    set(METIS_LIBRARY_DIR "${METIS_DIR}/lib" CACHE PATH "path to mpi library")
    set(METIS_LIBRARY_NAME "metis" CACHE STRING "name of metis library")

    link_directories(${METIS_LIBRARY_DIR})      #库目录
    include_directories(${METIS_HEADER_DIR})    #头文件目录
    link_libraries(${METIS_LIBRARY_NAME})       #库名称
endif()

#CGNS库设置
option(CGNS_ENABLE "ON for using CGNS library" ON)

if(CGNS_ENABLE)
    add_definitions(-DENABLE_CGNS)

    if(USE_PREBUILT_LIB)
        set(CGNS_DIR "${THIRD_LIB_DIR}/cgns")
        if( AUTO_CI_TEST )
            set(CGNS_DIR "C:/Program Files/cgns")
        endif()
    else()
        set(CGNS_DIR "/usr/local/cgns")
    endif()

    set(CGNS_HEADER_DIR "${CGNS_DIR}/include" CACHE PATH "path to cgns headers")
    set(CGNS_LIBRARY_DIR "${CGNS_DIR}/lib" CACHE PATH "path to cgns library")
    set(CGNS_LIBRARY_NAME "cgns" CACHE STRING "name of cgns library")

    include_directories(${CGNS_HEADER_DIR})  #头文件目录
    link_directories(${CGNS_LIBRARY_DIR})    #库目录
    link_libraries(${CGNS_LIBRARY_NAME})     #库名称
    
    #HDF5库
    #option(HDF5_ENABLE "ON for using HDF5 library" ON)
    set( HDF5_ENABLE TRUE CACHE BOOL "TRUE for using HDF5 library")
    set( HDF5_ENABLE_SZIP FALSE CACHE BOOL "TRUE for using HDF5 SZIP library")
    set( HDF5_ENABLE_ZLIB FALSE CACHE BOOL "TRUE for using HDF5 ZLIB library")
    if( AUTO_CI_TEST )
        #set( HDF5_ENABLE FALSE CACHE BOOL "FALSE for using HDF5 library" FORCE )
        set( HDF5_ENABLE TRUE CACHE BOOL "FALSE for using HDF5 library" FORCE )
    endif()

    if(HDF5_ENABLE)
        message(STATUS "The HDF5_ENABLE is ${HDF5_ENABLE}")
        add_definitions(-DENABLE_HDF5)

        if(USE_PREBUILT_LIB)
            set(HDF5_DIR_TMP "${THIRD_LIB_DIR}/hdf5")
            if( AUTO_CI_TEST )
                set(HDF5_DIR_TMP "C:/Program Files/HDF_Group/HDF5/1.10.5")
            endif()
            set(HDF5_LIBRARY_NAME "libhdf5" CACHE STRING "name of hdf5 library")
        else()
            set(HDF5_DIR_TMP "/usr/local/HDF_Group/HDF5/1.10.5")
            set(HDF5_LIBRARY_NAME "hdf5" CACHE STRING "name of hdf5 library")
        endif()

        set(HDF5_LIBRARY_SZIP_NAME "szip" CACHE STRING "name of szip library")
        set(HDF5_LIBRARY_ZLIB_NAME "zlib" CACHE STRING "name of zlib library")

        set(HDF5_DIR "${HDF5_DIR_TMP}" CACHE PATH "path to hdf5 directory" )

        set(HDF5_HEADER_DIR "${HDF5_DIR}/include" CACHE PATH "path to hdf5 headers" )
        set(HDF5_LIBRARY_DIR "${HDF5_DIR}/lib" CACHE PATH "path to hdf5 library" )

        #message( STATUS " The HDF5_DIR is ${HDF5_DIR} " )
        #message( STATUS " The HDF5_LIBRARY_NAME is ${HDF5_LIBRARY_NAME} " )
        #message( STATUS " The HDF5_HEADER_DIR is ${HDF5_HEADER_DIR} " )
        #message( STATUS " The HDF5_LIBRARY_DIR is ${HDF5_LIBRARY_DIR} " )

        include_directories(${HDF5_HEADER_DIR})
        link_directories(${HDF5_LIBRARY_DIR})
        link_libraries(${HDF5_LIBRARY_NAME})
        if ( HDF5_ENABLE_SZIP )
           link_libraries(${HDF5_LIBRARY_SZIP_NAME})
        endif()
        if ( HDF5_ENABLE_ZLIB )
           link_libraries(${HDF5_LIBRARY_ZLIB_NAME})
        endif()
    endif()
endif()

set( JsonPrjDir "${PRJ_ROOT_DIR}/${THIRD_PARTY_DIR}/jsoncpp" )
set( JsonHeaderDir ${JsonPrjDir} )
include_directories( ${JsonHeaderDir} ) 
set( JsonSRC "${JsonPrjDir}/jsoncpp.cpp" ) 

GetAllSubDir( dirList ${CMAKE_CURRENT_SOURCE_DIR})

message( STATUS "The dirList is ${dirList}" )

set_source_group( HEAD_LIST SRC_LIST )

list( APPEND SRC_LIST ${JsonSRC} ) 
#message(STATUS "SRC_LIST is ${SRC_LIST}.")

#. set executable name:
set(EXECUTABLE_TARGET_NAME "${CMAKE_PROJECT_NAME}")
message(STATUS "The executable name is ${EXECUTABLE_TARGET_NAME}.")

# adds flags to the compiler command line for sources in the current directory and below:
if( BIT64_FLAG )
    add_definitions(-DX64)
endif( BIT64_FLAG )

if( WIN32 )
    add_definitions( -DHX_PARALLEL -DHX_CGNS )
else()
    add_definitions( -DHX_PARALLEL -DHX_CGNS )
endif()

#. add an executable to the project using the specified source files and header files:
add_executable( ${EXECUTABLE_TARGET_NAME} ${SRC_LIST} ${HEAD_LIST} ) 

set(PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/)
set(PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/)

install(TARGETS ${EXECUTABLE_TARGET_NAME} RUNTIME DESTINATION bin)
install(DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system")
install(DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc")

#message(STATUS "The INCLUDE_DIRECTORIES name is ${INCLUDE_DIRECTORIES}.")
get_property( reqs TARGET ${EXECUTABLE_TARGET_NAME} PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "${EXECUTABLE_TARGET_NAME} requires ${reqs}")

message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
 