message ( STATUS "") 
message ( STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---" )
message ( STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}." )

set ( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code" )
set ( my_compile_def "" )
set ( my_libraries "" )
set ( my_include_dirs "" )

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
		list( APPEND my_compile_def WIN_GNU )
        set ( WIN_COMPILER mingw )
    else()
        set ( WIN_COMPILER ${VS_COMPILER} )
    endif()

    message ( STATUS "The current WIN_COMPILER is ${WIN_COMPILER}" )
endif()

if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
    message ( STATUS "WIN_COMPILER is ${WIN_COMPILER}" )
    #取消c语言中fopen等函数的警告
	list ( APPEND my_compile_def _CRT_SECURE_NO_WARNINGS )
endif()

message ( STATUS "The CMAKE_C_COMPILER is ${CMAKE_C_COMPILER}" )
message ( STATUS "The CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}" )
message ( STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}" )
message ( STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}" )

message ( STATUS "The CFLAGS is ${CFLAGS}" )
message ( STATUS "The CXXFLAGS is ${CXXFLAGS}" )
message ( STATUS "The LDFLAGS is ${LDFLAGS}" )

message ( STATUS "The CMAKE_C_FLAGS is ${CMAKE_C_FLAGS}" )
message ( STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}" )
message ( STATUS "The CMAKE_EXE_LINKER_FLAGS is ${CMAKE_EXE_LINKER_FLAGS}" )

if ( BIT64_FLAG )
    set ( ARCH "x64" )
else()
    set ( ARCH "x86" )
endif()

#mpi库设置
option ( MPI_ENABLE "ON for using MPI library" ON )

if ( MPI_ENABLE )
	list ( APPEND my_compile_def ENABLE_MPI )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
		set ( MPI_INCLUDE_DIRS "C:/Program Files (x86)/Microsoft SDKs/MPI/include" CACHE PATH "path to mpi headers" )
		set ( MPI_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${ARCH}/msmpi.lib" CACHE PATH "path to mpi library" )
    else()
		#set ( MPI_INCLUDE_DIRS "/usr/include/x86_64-linux-gnu/mpich" CACHE PATH "path to mpi headers" )
		#set ( MPI_LIBRARIES "/usr/lib/x86_64-linux-gnu/libmpich.so" CACHE PATH "path to mpi library" )
    set ( MPI_INCLUDE_DIRS "$ENV{MPI_HOME}/include" CACHE PATH "path to mpi headers" )
    set ( MPI_LIBRARIES "$ENV{MPI_HOME}/lib/libmpich.so" CACHE PATH "path to mpi library" )
    endif()
	
	message ( STATUS "The MPI_LIBRARIES is ${MPI_LIBRARIES}" )
	
	list ( APPEND my_libraries ${MPI_LIBRARIES} )
	list ( APPEND my_include_dirs ${MPI_INCLUDE_DIRS} )
endif()

#metis库设置
option ( METIS_ENABLE "ON for using METIS partitioning library" ON )

if ( METIS_ENABLE )
	list ( APPEND my_compile_def ENABLE_METIS )
    
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		set ( METIS_INCLUDE_DIRS "C:/METIS/include"  CACHE PATH "path to METIS headers" )
		set ( METIS_LIBRARIES "C:/METIS/lib/metis.lib" CACHE PATH "path to METIS library" )
    else()
		#set ( METIS_INCLUDE_DIRS "/usr/include" )
		#set ( METIS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libmetis.so" )
    set ( METIS_INCLUDE_DIRS "$ENV{METIS_HOME}/include" )
    set ( METIS_LIBRARIES "$ENV{METIS_HOME}/lib/libmetis.a" )
    endif()

	list ( APPEND my_libraries ${METIS_LIBRARIES} )
	list ( APPEND my_include_dirs ${METIS_INCLUDE_DIRS} )
endif()

#CGNS库设置
option ( CGNS_ENABLE "ON for using CGNS library" ON )

if ( CGNS_ENABLE )
	list ( APPEND my_compile_def ENABLE_CGNS )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
		set ( CGNS_INCLUDE_DIRS "C:/cgns/include" CACHE PATH "path to CGNS headers" )
		set ( CGNS_LIBRARIES "C:/cgns/lib/cgnsdll.lib" CACHE PATH "path to CGNS library" )
    else()
		#set ( CGNS_INCLUDE_DIRS "/usr/include" CACHE PATH "path to CGNS headers" )
		#set ( CGNS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libcgns.so" CACHE PATH "path to CGNS library" )
		set ( CGNS_INCLUDE_DIRS "$ENV{CGNS_HOME}/include" CACHE PATH "path to CGNS headers" )
		set ( CGNS_LIBRARIES "$ENV{CGNS_HOME}/lib/libcgns.so" CACHE PATH "path to CGNS library" )
    endif()
	message(STATUS "The CGNS_INCLUDE_DIRS is ${CGNS_INCLUDE_DIRS}")
	
	list ( APPEND my_libraries ${CGNS_LIBRARIES} )
	list ( APPEND my_include_dirs ${CGNS_INCLUDE_DIRS} )
endif()

set ( JsonPrjDir "${PRJ_ROOT_DIR}/${THIRD_PARTY_DIR}/jsoncpp" )
set ( JsonHeaderDir ${JsonPrjDir} )

list ( APPEND my_include_dirs ${JsonHeaderDir} )
set ( JsonSRC "${JsonPrjDir}/jsoncpp.cpp" ) 

ConstructSolutionDirTree( ${CMAKE_CURRENT_SOURCE_DIR} my_head_list my_src_list local_include_dirs )
list ( APPEND my_include_dirs ${local_include_dirs} )

list( APPEND my_src_list ${JsonSRC} ) 

if ( BIT64_FLAG )
    list( APPEND my_compile_def X64 )
endif( BIT64_FLAG )

if( WIN32 )
    list( APPEND my_compile_def HX_PARALLEL )
	list( APPEND my_compile_def HX_CGNS )
	list( APPEND my_compile_def USE_DLL )
else()
    list( APPEND my_compile_def HX_PARALLEL )
	list( APPEND my_compile_def HX_CGNS )
endif()

#. set executable name:
set ( OneFLOW_EXE "${CMAKE_PROJECT_NAME}" )
message ( STATUS "The executable name is ${OneFLOW_EXE}." )

#. add an executable to the project
add_executable ( ${OneFLOW_EXE} ${my_src_list} ${my_head_list} )

target_include_directories ( ${OneFLOW_EXE}
	PRIVATE 
        ${my_include_dirs}
)

target_link_libraries ( ${OneFLOW_EXE} 
    PRIVATE 
        ${my_libraries}
)

target_compile_definitions ( ${OneFLOW_EXE}
	PRIVATE
	   ${my_compile_def} 
)

set ( PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/ )
set ( PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/ )

install ( TARGETS ${OneFLOW_EXE} RUNTIME DESTINATION bin )
install ( DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system" )
install ( DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc" )

message ( STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---" )
