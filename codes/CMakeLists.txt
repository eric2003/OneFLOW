#. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code")
#set( HDF_FLAG FALSE CACHE BOOL "hdf or adf code")
set( HDF_FLAG TRUE CACHE BOOL "hdf or adf code")

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

#if (BIT64_FLAG)
#    set(BITCODE x64)
#else(BIT64_FLAG)
#    set(BITCODE x86)
#endif(BIT64_FLAG)

#if(WIN32)
#else()
#  set(CMAKE_CXX_COMPILER  mpicxx)
#  set(CMAKE_C_COMPILER      mpicc)
#endif(WIN32)

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

if(BIT64_FLAG)
    set(ARCH "x64")
else()
    set(ARCH "x86")
endif()

#使用已经编译好的库
if(USE_PREBUILT_LIB)
    set(THIRD_LIB_BASE_DIR "${CMAKE_SOURCE_DIR}/ThirdPartyLib")
    set(THIRD_LIB_DIR "${CMAKE_SOURCE_DIR}/ThirdPartyLib/${VS_COMPILER}/${ARCH}")
endif()

#mpi库设置(不需要区别vs版本)
option(MPI_ENABLE "ON for using MPI library" ON)

if(MPI_ENABLE)
    add_definitions(-DENABLE_MPI)
    
    if(USE_PREBUILT_LIB AND WIN32)
        set(MPI_HEADER_DIR "${THIRD_LIB_BASE_DIR}/msmpi/include" CACHE PATH "path to mpi headers")
        set(MPI_LIBRARY_DIR "${THIRD_LIB_BASE_DIR}/msmpi/lib/${ARCH}" CACHE PATH "path to mpi library")
    else()
        set(MPI_HEADER_DIR "/usr/local/mpich/include" CACHE PATH "path to mpi headers")
        set(MPI_LIBRARY_DIR "/usr/local/mpich/lib" CACHE PATH "path to mpi library")
    endif()

    if( AUTO_CI_TEST )
        # determine if MPI_DIR has been set
        #-DMPI_DIR=$TRAVIS_BUILD_DIR/openmpi
        if (MPI_DIR)
            message (STATUS "Looking for MPI here: ${MPI_DIR}")
            set (CMAKE_PREFIX_PATH ${MPI_DIR})
        endif ()

        find_package (MPI REQUIRED)
        include_directories (SYSTEM "${MPI_INCLUDE_PATH}")

        message( STATUS "MPI_DIR is ${MPI_DIR}" )
        message( STATUS "MPI_INCLUDE_PATH is ${MPI_INCLUDE_PATH}" )
        message( STATUS "MPI_C_LIBRARIES is ${MPI_C_LIBRARIES}" )

        message( STATUS "MPI_HEADER_DIR is ${MPI_HEADER_DIR}" )
        message( STATUS "MPI_LIBRARY_DIR is ${MPI_LIBRARY_DIR}" )
    endif ()

    if(WIN32)
        set(MPI_LIBRARY_NAME "msmpi" CACHE STRING "name of mpi library")
    else()
        set(MPI_LIBRARY_NAME "mpi" CACHE STRING "name of mpi library")
    endif()

    include_directories(${MPI_HEADER_DIR})  #头文件目录
    link_directories(${MPI_LIBRARY_DIR})    #库目录
    link_libraries(${MPI_LIBRARY_NAME})     #库名称
endif()

#metis库设置
option(METIS_ENABLE "ON for using METIS partitioning library" ON)

if(METIS_ENABLE)
    add_definitions(-DENABLE_METIS)
    
    if(USE_PREBUILT_LIB)
        set(METIS_HEADER_DIR "${THIRD_LIB_DIR}/metis/include" CACHE PATH "path to metis headers")
        set(METIS_LIBRARY_DIR "${THIRD_LIB_DIR}/metis/lib" CACHE PATH "path to mpi library")
    else()
        set(METIS_HEADER_DIR "/usr/local/metis/include" CACHE PATH "path to metis headers")
        set(METIS_LIBRARY_DIR "/usr/local/metis/lib" CACHE PATH "path to mpi library")
    endif()

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
        set(CGNS_HEADER_DIR "${THIRD_LIB_DIR}/cgns/include" CACHE PATH "path to cgns headers")
        set(CGNS_LIBRARY_DIR "${THIRD_LIB_DIR}/cgns/lib" CACHE PATH "path to cgns library")
    else()
        set(CGNS_HEADER_DIR "/usr/local/cgns/include" CACHE PATH "path to cgns headers")
        set(CGNS_LIBRARY_DIR "/usr/local/cgns/lib" CACHE PATH "path to cgns library")
    endif()
    
    set(CGNS_LIBRARY_NAME "cgns" CACHE STRING "name of cgns library")

    include_directories(${CGNS_HEADER_DIR})  #头文件目录
    link_directories(${CGNS_LIBRARY_DIR})    #库目录
    link_libraries(${CGNS_LIBRARY_NAME})     #库名称
    
    #HDF5库
    option(HDF5_ENABLE "ON for using HDF5 library" ON)

    if(HDF5_ENABLE)
        add_definitions(-DENABLE_HDF5)

        if(USE_PREBUILT_LIB)
            set(HDF5_HEADER_DIR "${THIRD_LIB_DIR}/hdf5/include" CACHE PATH "path to hdf5 headers")
            set(HDF5_LIBRARY_DIR "${THIRD_LIB_DIR}/hdf5/lib" CACHE PATH "path to hdf5 library")
            set(HDF5_LIBRARY_NAME "libhdf5" CACHE STRING "name of hdf5 library")
        else()
            set(HDF5_HEADER_DIR "/usr/local/HDF_Group/HDF5/1.10.5/include" CACHE PATH "path to hdf5 headers")
            set(HDF5_LIBRARY_DIR "/usr/local/HDF_Group/HDF5/1.10.5/lib" CACHE PATH "path to hdf5 library")
            set(HDF5_LIBRARY_NAME "hdf5" CACHE STRING "name of hdf5 library")
        endif()

        include_directories(${HDF5_HEADER_DIR})
        link_directories(${HDF5_LIBRARY_DIR})
        link_libraries(${HDF5_LIBRARY_NAME})
    endif()
endif()

#if(WIN32)
#    link_directories("${PROJECT_SOURCE_DIR}/lib/windows/${VS_COMPILER}/${BITCODE}/")
#else(WIN32)
#    link_directories("${PROJECT_SOURCE_DIR}/lib/linux/${BITCODE}/")
#endif(WIN32)

function(GetALLSUBDIR result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist} PARENT_SCOPE)
endfunction()

GetALLSUBDIR( dirList ${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "The executable name is ${dirList}" )

foreach(dir ${dirList})
    message( ${dir} )
    file(GLOB headerFile "${dir}/include/*.h" )
    file(GLOB srcFile "${dir}/src/*.cpp" )
    file(GLOB hppFile "${dir}/src/*.hpp" )
    message( STATUS ${srcFile} )
    message( STATUS ${headerFile} )
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/${dir}/include")
    if( NOT( hppFile STREQUAL "" ) )
        include_directories("${CMAKE_CURRENT_SOURCE_DIR}/${dir}/src")
    endif()

    source_group("${dir}\\include" FILES ${headerFile})
    source_group("${dir}\\src" FILES ${srcFile})
    list(APPEND HEAD_LIST ${headerFile})
    list(APPEND SRC_LIST ${srcFile})
endforeach()


#. set executable name:
set(EXECUTABLE_TARGET_NAME "${CMAKE_PROJECT_NAME}")
message(STATUS "The executable name is ${EXECUTABLE_TARGET_NAME}.")

# adds flags to the compiler command line for sources in the current directory and below:
if ( BIT64_FLAG )
     add_definitions(-DX64)
endif( BIT64_FLAG )

if(WIN32)
   add_definitions(-DHX_PARALLEL -DHX_CGNS)
else(WIN32)
   add_definitions(-DHX_PARALLEL -DHX_CGNS)
endif(WIN32)

#. add an executable to the project using the specified source files and header files:
add_executable(${EXECUTABLE_TARGET_NAME} ${SRC_LIST} ${HEAD_LIST}) 

#if(WIN32)
#  if(HDF_FLAG)
#    target_link_libraries(${EXECUTABLE_TARGET_NAME} msmpi cgns libhdf5 metis)
#  else(HDF_FLAG)
#    target_link_libraries(${EXECUTABLE_TARGET_NAME} msmpi cgns metis)
#  endif(HDF_FLAG)
#else(WIN32)
#    if(HDF_FLAG)
#        target_link_libraries(${EXECUTABLE_TARGET_NAME} libmetis.a libcgns.a libhdf5.a )
#    else(HDF_FLAG)
#         target_link_libraries(${EXECUTABLE_TARGET_NAME} libmetis.a libcgns.a )
#    endif(HDF_FLAG)
#endif(WIN32)   
                   
#3. exit:          
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
 