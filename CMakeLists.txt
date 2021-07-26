cmake_minimum_required ( VERSION 3.15 )

message(STATUS "The CMAKE_VERSION is ${CMAKE_VERSION}.")
include( "${CMAKE_CURRENT_SOURCE_DIR}/project-meta-info.in" )
		
project ( OneFLOW
          VERSION ${project_version}
          DESCRIPTION ${project_description}
          HOMEPAGE_URL ${project_homepage}
          LANGUAGES CXX
	     )

message ( STATUS "CMAKE_PROJECT_NAME    = ${CMAKE_PROJECT_NAME}" )
message ( STATUS "PROJECT_VERSION_MAJOR = ${PROJECT_VERSION_MAJOR}" )
message ( STATUS "PROJECT_VERSION_MINOR = ${PROJECT_VERSION_MINOR}" )
message ( STATUS "PROJECT_VERSION_PATCH = ${PROJECT_VERSION_PATCH}" )
message ( STATUS "PROJECT_VERSION_TWEAK = ${PROJECT_VERSION_TWEAK}" )
message ( STATUS "PROJECT_HOMEPAGE_URL  = ${PROJECT_HOMEPAGE_URL}")
message ( STATUS "CMAKE_PROJECT_DESCRIPTION = ${CMAKE_PROJECT_DESCRIPTION}" )

#set working directory:
set ( MY_WORKING_DIR OneFLOW     )
set ( THIRD_PARTY_DIR ThirdParty )

set ( PRJ_ROOT_DIR ${CMAKE_SOURCE_DIR} )
message ( STATUS "PRJ_ROOT_DIR is ${PRJ_ROOT_DIR}." )


set ( OneFLOW_CMAKE_DIR "${CMAKE_SOURCE_DIR}/cmake" )
list ( INSERT CMAKE_MODULE_PATH 0 "${OneFLOW_CMAKE_DIR}" )

#output message:
message ( STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}." )
message ( STATUS "The OS is ${CMAKE_SYSTEM_NAME}, version ${CMAKE_SYSTEM_VERSION}." )

message ( STATUS "The CMAKE_SOURCE_DIR is ${CMAKE_SOURCE_DIR}." )
message ( STATUS "The CMAKE_BINARY_DIR is ${CMAKE_BINARY_DIR}." ) 

message ( STATUS "The working directory is ${CMAKE_BINARY_DIR}${MY_WORKING_DIR}." )

message ( STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---" ) 
message ( STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}." )

message ( STATUS "The C compiler is ${CMAKE_C_COMPILER_ID}" )
message ( STATUS "The CXX compiler is ${CMAKE_CXX_COMPILER_ID}" )
message ( STATUS "The CMAKE_GENERATOR is ${CMAKE_GENERATOR}" )
message ( STATUS "The CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}" )

message ( STATUS "The CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}" )
 
set ( CMAKE_CXX_STANDARD 11 )
set ( CMAKE_CXX_STANDARD_REQUIRED True )
set ( CMAKE_CXX_EXTENSIONS OFF )

include ( UtilFunction   )
include ( oflDirectories )

# Set up our directory structure for output libraries and binaries
set ( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" )
set ( CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" )
set ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" )

message ( STATUS "The CMAKE_BINARY_DIR is ${CMAKE_BINARY_DIR}" )
message ( STATUS "The CMAKE_INSTALL_BINDIR is ${CMAKE_INSTALL_BINDIR}" )
message ( STATUS "The CMAKE_INSTALL_LIBDIR is ${CMAKE_INSTALL_LIBDIR}" )
message ( STATUS "The CMAKE_RUNTIME_OUTPUT_DIRECTORY is ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" )
message ( STATUS "The CMAKE_LIBRARY_OUTPUT_DIRECTORY is ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" )
message ( STATUS "The CMAKE_ARCHIVE_OUTPUT_DIRECTORY is ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}" )

#look for CMakeLists.txt in subdirectory <codes> and build an executable target:
add_subdirectory ( ThirdParty ${THIRD_PARTY_DIR} )
add_subdirectory ( codes      ${MY_WORKING_DIR}  )

#exit:
message ( STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---" )
message ( STATUS "Finished build project ${PROJECT_NAME}!" )
