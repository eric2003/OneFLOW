message ( STATUS "") 
message ( STATUS "--- Enter subdirectory, Start to call ${CMAKE_CURRENT_LIST_FILE} ---" )
message ( STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}." )

set ( my_compile_def )
set ( my_libraries )
set ( my_include_dirs )
set ( my_propertie_name )
set ( my_propertie_value )
set ( my_compile_features )

list ( APPEND my_compile_features cxx_std_14 )

set ( BIT64_FLAG TRUE CACHE BOOL "64bit or 32bit code" )

include ( CheckLanguage )
check_language ( CUDA )

if ( CMAKE_CUDA_COMPILER )
	enable_language ( CUDA )
	list( APPEND my_compile_def ENABLE_CUDA )
    #set ( PRJ_CUDA_ARCHITECTURES "35;50;52;72;75" )
	#set ( PRJ_CUDA_SEPARABLE_COMPILATION ON )
	list( APPEND my_include_dirs ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} )
	message ( STATUS "The CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES is ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}" )
	list ( APPEND my_propertie_name CUDA_ARCHITECTURES )
	list ( APPEND my_propertie_value "35;50;52;72;75" )
	list ( APPEND my_propertie_name CUDA_SEPARABLE_COMPILATION )
	list ( APPEND my_propertie_value ON )
	list ( APPEND my_compile_features cuda_std_14 )
else()
	message ( STATUS "No CUDA support" )
endif()

message ( STATUS "The CMAKE_CUDA_COMPILER is ${CMAKE_CUDA_COMPILER}" )


if ( CMAKE_SYSTEM_NAME MATCHES "Linux" )
elseif ( CMAKE_SYSTEM_NAME MATCHES "Windows" )
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
        #set ( MPI_INCLUDE_DIRS "$ENV{MPI_HOME}/include" CACHE PATH "path to mpi headers" )
        #set ( MPI_LIBRARIES "$ENV{MPI_HOME}/lib/libmpich.so" CACHE PATH "path to mpi library" )
        set ( MPI_INCLUDE_DIRS "$ENV{MPI_HOME_INC}" CACHE PATH "path to mpi headers" )
        set ( MPI_LIBRARIES "$ENV{MPI_HOME_LIB}" CACHE PATH "path to mpi library" )
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
		#set ( METIS_INCLUDE_DIRS "$ENV{METIS_HOME}/include" )
		#set ( METIS_LIBRARIES "$ENV{METIS_HOME}/lib/libmetis.a" )
		set ( METIS_INCLUDE_DIRS $ENV{METIS_HOME_INC} )
		set ( METIS_LIBRARIES $ENV{METIS_HOME_LIB} )
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
		#set ( CGNS_INCLUDE_DIRS "$ENV{CGNS_HOME}/include" CACHE PATH "path to CGNS headers" )
		#set ( CGNS_LIBRARIES "$ENV{CGNS_HOME}/lib/libcgns.so" CACHE PATH "path to CGNS library" )
		set ( CGNS_INCLUDE_DIRS $ENV{CGNS_HOME_INC} CACHE PATH "path to CGNS headers" )
		set ( CGNS_LIBRARIES $ENV{CGNS_HOME_LIB} CACHE PATH "path to CGNS library" )		
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
set ( OneFLOW_EXE "${PROJECT_NAME}" )
message ( STATUS "The executable name is ${PROJECT_NAME}." )

#. add an executable to the project
add_executable ( ${PROJECT_NAME} ${my_src_list} ${my_head_list} )

target_include_directories ( ${PROJECT_NAME}
	PRIVATE 
        ${my_include_dirs}
)

target_link_libraries ( ${PROJECT_NAME} 
    PRIVATE 
        ${my_libraries}
)

target_compile_features ( ${PROJECT_NAME} 
	PUBLIC 
		${my_compile_features}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${my_compile_def} 
)

list ( LENGTH my_propertie_name len )
math ( EXPR my_len "${len} - 1")

foreach ( val RANGE ${my_len} )
	list ( GET my_propertie_name  ${val} name )
	list ( GET my_propertie_value ${val} value )
	set_target_properties ( ${PROJECT_NAME}
		PROPERTIES
			${name} ${value}
	)
endforeach ()

#set_target_properties (	${PROJECT_NAME}
#	PROPERTIES
#		CUDA_ARCHITECTURES ${PRJ_CUDA_ARCHITECTURES}
#		CUDA_SEPARABLE_COMPILATION ${PRJ_CUDA_SEPARABLE_COMPILATION}
#)

set ( PRJ_SYSTEM_DIR ${PRJ_ROOT_DIR}/system/ )
set ( PRJ_DOCUMENT_DIR ${PRJ_ROOT_DIR}/doc/ )

install ( TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin )
install ( DIRECTORY "${PRJ_SYSTEM_DIR}" DESTINATION "bin/system" )
install ( DIRECTORY "${PRJ_DOCUMENT_DIR}" DESTINATION "doc" )

message ( STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---" )
