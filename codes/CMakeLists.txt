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

if (BIT64_FLAG)
    set(BITCODE x64)
else(BIT64_FLAG)
    set(BITCODE x86)
endif(BIT64_FLAG)

if(WIN32)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER      mpicc)
endif(WIN32)

message(STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}")
message(STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message(STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}")
message(STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

if(WIN32)
    link_directories("${PROJECT_SOURCE_DIR}/lib/windows/${VS_COMPILER}/${BITCODE}/")
else(WIN32)
    link_directories("${PROJECT_SOURCE_DIR}/lib/linux/${BITCODE}/")
endif(WIN32)

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

if(WIN32)
  if(HDF_FLAG)
    target_link_libraries(${EXECUTABLE_TARGET_NAME} msmpi cgns libhdf5 metis)
  else(HDF_FLAG)
    target_link_libraries(${EXECUTABLE_TARGET_NAME} msmpi cgns metis)
  endif(HDF_FLAG)
else(WIN32)
    if(HDF_FLAG)
        target_link_libraries(${EXECUTABLE_TARGET_NAME} libmetis.a libcgns.a libhdf5.a )
    else(HDF_FLAG)
         target_link_libraries(${EXECUTABLE_TARGET_NAME} libmetis.a libcgns.a )
    endif(HDF_FLAG)
endif(WIN32)   
                   
#3. exit:          
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
 