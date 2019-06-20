#. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code")
#set( HDF_FLAG FALSE CACHE BOOL "hdf or adf code")
set( HDF_FLAG TRUE CACHE BOOL "hdf or adf code")

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

if(WIN32)
message(STATUS "VS_COMPILER is ${VS_COMPILER}")
endif(WIN32)

#if (BIT64_FLAG)
#    set(BITCODE x64)
#else(BIT64_FLAG)
#    set(BITCODE x86)
#endif(BIT64_FLAG)

if(WIN32)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER      mpicc)
endif(WIN32)

message(STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}")
message(STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message(STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message(STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

#metis库设置
option(METIS_ENABLE "ON for METIS partitioning" ON)

if(METIS_ENABLE)
	add_definitions(-DENABLE_METIS)
	
	set(METIS_HEADER_DIR "/path/to/metis/headers" CACHE PATH "path to metis headers")
	include_directories(${METIS_HEADER_DIR})	#头文件目录

	set(METIS_LIBRARY_DIR "/path/to/metis/library" CACHE PATH "path to mpi library")
	link_directories(${METIS_LIBRARY_DIR})		#库目录
	
	set(METIS_LIBRARY_NAME "metis" CACHE STRING "name of metis library")
	link_libraries(${METIS_LIBRARY_NAME})		#库名称
endif()

#mpi库设置
option(MPI_ENABLE "ON for MPI parallel" ON)

if(MPI_ENABLE)
	add_definitions(-DENABLE_MPI)
	
	set(MPI_HEADER_DIR "/path/to/mpi/headers" CACHE PATH "path to mpi headers")
	include_directories(${MPI_HEADER_DIR})	#头文件目录

	set(MPI_LIBRARY_DIR "/path/to/mpi/library" CACHE PATH "path to mpi library")
	link_directories(${MPI_LIBRARY_DIR})		#库目录
	
	set(MPI_LIBRARY_NAME "mpi" CACHE STRING "name of mpi library")
	link_libraries(${MPI_LIBRARY_NAME})			#库名称
endif()


#CGNS库设置
option(CGNS_ENABLE "ON for CGNS" ON)

if(CGNS_ENABLE)
	add_definitions(-DENABLE_CGNS)
	
	set(CGNS_HEADER_DIR "/path/to/cgns/headers" CACHE PATH "path to cgns headers")
	include_directories(${CGNS_HEADER_DIR})	#头文件目录

	set(CGNS_LIBRARY_DIR "/path/to/cgns/library" CACHE PATH "path to cgns library")
	link_directories(${CGNS_LIBRARY_DIR})		#库目录
	
	set(CGNS_LIBRARY_NAME "cgns" CACHE STRING "name of cgns library")
	link_libraries(${CGNS_LIBRARY_NAME})		#库名称
	
	#HDF5库
	option(HDF5_ENABLE "ON for HDF5" ON)

	if(HDF5_ENABLE)
		add_definitions(-DENABLE_HDF5)
		
		#头文件目录
		set(HDF5_HEADER_DIR "/path/to/hdf5/headers" CACHE PATH "path to hdf5 headers")
		include_directories(${HDF5_HEADER_DIR})	
		#库目录
		set(HDF5_LIBRARY_DIR "/path/to/hdf5/library" CACHE PATH "path to hdf5 library")
		link_directories(${HDF5_LIBRARY_DIR})		
		#库名称
		set(HDF5_LIBRARY_NAME "libhdf5" CACHE STRING "name of hdf5 library")
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
 