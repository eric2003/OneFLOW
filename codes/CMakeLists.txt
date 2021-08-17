message ( STATUS "") 
message ( STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---" )
message ( STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}." )

set ( PRJ_COMPILE_DEF )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )
set ( PRJ_PROPERTY_NAME )
set ( PRJ_PROPERTY_VALUE )
set ( PRJ_COMPILE_FEATURES )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_14 )

set ( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code" )

include ( CheckLanguage )
check_language ( CUDA )

if ( CMAKE_CUDA_COMPILER )
	enable_language ( CUDA )
	list ( APPEND PRJ_COMPILE_DEF ENABLE_CUDA )
	list ( APPEND PRJ_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} )
	list ( APPEND PRJ_PROPERTY_NAME CUDA_ARCHITECTURES )
	list ( APPEND PRJ_PROPERTY_VALUE "35;50;52;72;75" )
	list ( APPEND PRJ_PROPERTY_NAME CUDA_SEPARABLE_COMPILATION )
	list ( APPEND PRJ_PROPERTY_VALUE ON )
	list ( APPEND PRJ_COMPILE_FEATURES cuda_std_14 )
else()
	message ( STATUS "No CUDA support" )
endif()

message ( STATUS "The CMAKE_CUDA_COMPILER is ${CMAKE_CUDA_COMPILER}" )

if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
    #Cancel the warning of fopen in C language
	list ( APPEND PRJ_COMPILE_DEF _CRT_SECURE_NO_WARNINGS )
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

#MPI library settings
option ( MPI_ENABLE "ON for using MPI library" ON )

if ( MPI_ENABLE )
	list ( APPEND PRJ_COMPILE_DEF ENABLE_MPI )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
		set ( MPI_INCLUDE_DIRS "C:/Program Files (x86)/Microsoft SDKs/MPI/include" CACHE PATH "path to mpi headers" )
		set ( MPI_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/${ARCH}/msmpi.lib" CACHE PATH "path to mpi library" )
    else()
		#set ( MPI_INCLUDE_DIRS "/usr/include/x86_64-linux-gnu/mpich" CACHE PATH "path to mpi headers" )
		#set ( MPI_LIBRARIES "/usr/lib/x86_64-linux-gnu/libmpich.so" CACHE PATH "path to mpi library" )
        #set ( MPI_INCLUDE_DIRS "$ENV{MPI_HOME}/include" CACHE PATH "path to mpi headers" )
        #set ( MPI_LIBRARIES "$ENV{MPI_HOME}/lib/libmpich.so" CACHE PATH "path to mpi library" )
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
		set ( METIS_INCLUDE_DIRS "C:/METIS/include"  CACHE PATH "path to METIS headers" )
		set ( METIS_LIBRARIES "C:/METIS/lib/metis.lib" CACHE PATH "path to METIS library" )
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
		set ( CGNS_INCLUDE_DIRS "C:/cgns/include" CACHE PATH "path to CGNS headers" )
		set ( CGNS_LIBRARIES "C:/cgns/lib/cgnsdll.lib" CACHE PATH "path to CGNS library" )
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

list ( APPEND PRJ_INCLUDE_DIRS ${JsonHeaderDir} )
set ( JsonSRC "${JsonPrjDir}/jsoncpp.cpp" ) 

ConstructSolutionDirTree( ${CMAKE_CURRENT_SOURCE_DIR} PRJ_HEADER_LIST PRJ_SRC_LIST PRJ_LOCAL_INCLUDE_DIRS )
list ( APPEND PRJ_INCLUDE_DIRS ${PRJ_LOCAL_INCLUDE_DIRS} )

list ( APPEND PRJ_SRC_LIST ${JsonSRC} ) 

if ( BIT64_FLAG )
    list( APPEND PRJ_COMPILE_DEF X64 )
endif( BIT64_FLAG )

if( WIN32 )
    list( APPEND PRJ_COMPILE_DEF HX_PARALLEL )
	list( APPEND PRJ_COMPILE_DEF HX_CGNS )
	list( APPEND PRJ_COMPILE_DEF USE_DLL )
else()
    list( APPEND PRJ_COMPILE_DEF HX_PARALLEL )
	list( APPEND PRJ_COMPILE_DEF HX_CGNS )
endif()

#. set executable name:
set ( OneFLOW_EXE "${PROJECT_NAME}" )
message ( STATUS "The executable name is ${PROJECT_NAME}." )

#. add an executable to the project
add_executable ( ${PROJECT_NAME} ${PRJ_SRC_LIST} ${PRJ_HEADER_LIST} )

target_include_directories ( ${PROJECT_NAME}
	PRIVATE 
        ${PRJ_INCLUDE_DIRS}
)

target_link_libraries ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
	PUBLIC 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_COMPILE_DEF} 
)

list ( LENGTH PRJ_PROPERTY_NAME len )
math ( EXPR my_len "${len} - 1")

foreach ( val RANGE ${my_len} )
	list ( GET PRJ_PROPERTY_NAME  ${val} name )
	list ( GET PRJ_PROPERTY_VALUE ${val} value )
	set_target_properties ( ${PROJECT_NAME}
		PROPERTIES
			${name} ${value}
	)
endforeach ()

set ( PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/ )
set ( PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/ )

install ( TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin )
install ( DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system" )
install ( DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc" )

message ( STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---" )
