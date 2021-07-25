#. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code")

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    if ( CMAKE_GENERATOR MATCHES "Visual Studio 16 2019" )
        set( VS_COMPILER vs2019 )
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 15 2017" )
        set( VS_COMPILER vs2017 )
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 14 2015" )
        set( VS_COMPILER vs2015 )
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 12 2013" )
        set( VS_COMPILER vs2013 )
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 11 2012" )
        set( VS_COMPILER vs2012 )
    elseif (CMAKE_GENERATOR MATCHES "Visual Studio 10 2010" )
        set( VS_COMPILER vs2010 )
    else()
        set( VS_COMPILER vs2019 )
    endif()

    if ( CMAKE_COMPILER_IS_GNUCXX )
        add_definitions( -DWIN_GNU )
    endif()

    if ( CMAKE_COMPILER_IS_GNUCXX )
        add_definitions( -DWIN_GNU )
        set( WIN_COMPILER mingw )
    else()
        set( WIN_COMPILER ${VS_COMPILER} )
    endif()

    message(STATUS "The current WIN_COMPILER is ${WIN_COMPILER}")
endif()

if( CMAKE_SYSTEM_NAME MATCHES "Linux" )
elseif( CMAKE_SYSTEM_NAME MATCHES "Windows" )
    message(STATUS "WIN_COMPILER is ${WIN_COMPILER}")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)  #取消c语言中fopen等函数的警告
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

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if ( CMAKE_COMPILER_IS_GNUCXX )
        #Set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-Bstatic, -Wl,-Bdynamic")
    else()
        #Set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl")
    endif()
else()
    #Set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -ldl")
endif()
#message(STATUS "The CMAKE_EXE_LINKER_FLAGS is ${CMAKE_EXE_LINKER_FLAGS}")

if( BIT64_FLAG )
    set(ARCH "x64")
else()
    set(ARCH "x86")
endif()

#mpi库设置
option(MPI_ENABLE "ON for using MPI library" ON)

if ( MPI_ENABLE )
    add_definitions(-DENABLE_MPI)
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		set( MPI_INCLUDE_DIRS "C:/Program Files (x86)/Microsoft SDKs/MPI/include" CACHE PATH "path to mpi headers" )
		set( MPI_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${ARCH}/msmpi.lib" CACHE PATH "path to mpi library" )
    else()
		set( MPI_INCLUDE_DIRS "/usr/include/x86_64-linux-gnu/mpich" CACHE PATH "path to mpi headers" )
		set( MPI_LIBRARIES "/usr/lib/x86_64-linux-gnu/libmpich.so" CACHE PATH "path to mpi library" )
    endif()

	include_directories ("${MPI_INCLUDE_DIRS}")
	message(STATUS "The MPI_LIBRARIES is ${MPI_LIBRARIES}")
endif()

#metis库设置
option(METIS_ENABLE "ON for using METIS partitioning library" ON)

if( METIS_ENABLE )
    add_definitions( -DENABLE_METIS )
    
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		set( METIS_INCLUDE_DIRS "C:/METIS/include"  CACHE PATH "path to METIS headers" )
		set( METIS_LIBRARIES "C:/METIS/lib/metis.lib" CACHE PATH "path to METIS library" )
    else()
		set( METIS_INCLUDE_DIRS "/usr/include" )
		set( METIS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libmetis.so" )
    endif()

    include_directories(${METIS_INCLUDE_DIRS})
endif()

#CGNS库设置
option( CGNS_ENABLE "ON for using CGNS library" ON )

if( CGNS_ENABLE )
    add_definitions( -DENABLE_CGNS )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
		set( CGNS_INCLUDE_DIRS "C:/cgns/include" CACHE PATH "path to CGNS headers" )
		set( CGNS_LIBRARY_DIRS "C:/cgns/lib/cgnsdll.lib" CACHE PATH "path to CGNS library" )
    else()
		set( CGNS_INCLUDE_DIRS "/usr/include" CACHE PATH "path to CGNS headers" )
		set( CGNS_LIBRARY_DIRS "/usr/lib/x86_64-linux-gnu/libcgns.so" CACHE PATH "path to CGNS library" )
    endif()
	
    include_directories( ${CGNS_INCLUDE_DIRS} )
	message(STATUS "The CGNS_INCLUDE_DIRS is ${CGNS_INCLUDE_DIRS}")
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
    add_definitions( -DHX_PARALLEL -DHX_CGNS -DUSE_DLL )
else()
    add_definitions( -DHX_PARALLEL -DHX_CGNS )
endif()

#. add an executable to the project using the specified source files and header files:
add_executable( ${EXECUTABLE_TARGET_NAME} ${SRC_LIST} ${HEAD_LIST} ) 
target_link_libraries(${EXECUTABLE_TARGET_NAME} ${CGNS_LIBRARY_DIRS} ${MPI_LIBRARIES} ${METIS_LIBRARIES} )

set(PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/)
set(PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/)

install(TARGETS ${EXECUTABLE_TARGET_NAME} RUNTIME DESTINATION bin)
install(DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system")
install(DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc")

message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
