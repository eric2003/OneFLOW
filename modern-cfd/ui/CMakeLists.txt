cmake_minimum_required ( VERSION 3.20 )

#set( VCPKG_ROOT "c:/dev/vcpkg" )
#set( CMAKE_TOOLCHAIN_FILE "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" )

#project( OneFLOW-UI )
set ( GUI_NAME "OneFLOW-UI" )

#message( STATUS "VCPKG_ROOT=${VCPKG_ROOT}")
#message( STATUS "CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}")
#message( STATUS "PROJECT_NAME=${PROJECT_NAME}")

set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

find_package ( Qt6 COMPONENTS Widgets )
find_package ( Qt6 REQUIRED COMPONENTS Sql )

#list ( APPEND PRJ_LIBRARIES Qt6::Core )
list ( APPEND PRJ_LIBRARIES Qt6::Sql )
list ( APPEND PRJ_LIBRARIES Qt6::Widgets )

#list ( APPEND PRJ_INCLUDE_DIRS "c:\Program Files\MySQL\MySQL Server 8.0\include" )
list ( APPEND PRJ_LIBRARIES "c:/Program Files/MySQL/MySQL Server 8.0/lib/libmysql.lib" )


get_directory_property( my_system_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} BUILDSYSTEM_TARGETS )
get_directory_property( my_import_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} IMPORTED_TARGETS )

message( STATUS "my_system_targets=${my_system_targets}" )
message( STATUS "my_import_targets=${my_import_targets}" )

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${GUI_NAME} )
endif()

add_executable ( ${GUI_NAME}
    main.cpp
    mainwindow.cpp
    mainwindow.ui
	MyDataBase.cpp
	MyDataBase.h
)

set_property( TARGET ${GUI_NAME} PROPERTY AUTOUIC ON )
set_property( TARGET ${GUI_NAME} PROPERTY AUTOMOC ON )
set_property( TARGET ${GUI_NAME} PROPERTY AUTORCC ON )

target_include_directories ( ${GUI_NAME}
    PRIVATE
        ${PRJ_INCLUDE_DIRS}
)

target_link_libraries( ${GUI_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)

target_compile_features ( ${GUI_NAME} 
    PRIVATE 
        ${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${GUI_NAME}
    PRIVATE
       ${PRJ_COMPILE_DEFINITIONS} 
)

get_target_property( my_location Qt6::windeployqt IMPORTED_LOCATION )
get_target_property( my_location_exe Qt6::windeployqt IMPORTED_LOCATION_<CONFIG> )
get_target_property( my_location_debug Qt6::windeployqt IMPORTED_LOCATION_DEBUG )
get_target_property( my_location_release Qt6::windeployqt IMPORTED_LOCATION_RELEASE )

message ( STATUS "my_location = ${my_location}" )
message ( STATUS "my_location_exe = ${my_location_exe}" )
message ( STATUS "my_location_debug = ${my_location_debug}" )
message ( STATUS "my_location_release = ${my_location_release}" )
message ( STATUS "OUTPUT_BIN_DIR = ${OUTPUT_BIN_DIR}" )
message ( STATUS "OUTPUT_BIN_DIR/GUI_NAME = ${OUTPUT_BIN_DIR}/${GUI_NAME}" )

#add_custom_command ( TARGET ${GUI_NAME}
#                    POST_BUILD
#					COMMAND echo Qt6::windeployqt=$<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}>
#					COMMAND echo "Qt6::windeployqt=$<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}>"
#					COMMAND echo TARGET_FILE:${GUI_NAME}=$<TARGET_FILE:${GUI_NAME}>
#					COMMAND echo TARGET_FILE:${GUI_NAME}="$<TARGET_FILE:${GUI_NAME}>"
#					COMMAND $<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}> "$<TARGET_FILE:${GUI_NAME}>"
#					COMMAND $<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}> $<TARGET_FILE:${GUI_NAME}>
#                    COMMENT "Run Qt6::windeployqt" )

add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
					COMMAND echo Qt6::windeployqt=$<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}>
					#COMMAND C:/dev/vcpkg/installed/x64-windows/tools/Qt6/bin/windeployqt.exe $<TARGET_FILE:${GUI_NAME}>
					COMMAND $<$<CONFIG:Debug>:${my_location_debug}>$<$<CONFIG:Release>:${my_location_release}> $<TARGET_FILE:${GUI_NAME}>
                    COMMENT "Run MyQt6::windeployqt" )


message ( STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}" )
message ( STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}" )
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}" )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}" )

set ( PRJ_BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}" )
set ( PRJ_BIN_DIR1 "${CMAKE_CURRENT_BINARY_DIR}" )
message ( STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY_Debug = ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_Debug}" )
message ( STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY_Release = ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_Release}" )

if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    #set ( PRJ_BIN_DIR "${PRJ_BIN_DIR}/$<CONFIG>" )
endif()
message ( STATUS "PRJ_BIN_DIR = ${PRJ_BIN_DIR}" )

add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo "RUN POST-BUILD CUSTOM COMMAND"
                    COMMENT "Run custom command" )
					
add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo "RUN POST-BUILD CUSTOM COMMAND111111111111"
                    COMMENT "Run custom command" )

add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo CONFIG=$<CONFIG>
                    COMMENT "Run custom command" )
					
add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo my_location_debug=${my_location_debug}
                    COMMENT "Run custom command" )
					
add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo my_location_exe=${my_location_exe}
					COMMAND echo my_location_config=$<$<CONFIG:Debug>:${my_location_debug}> $<$<CONFIG:Release>:${my_location_release}>
                    COMMENT "Run custom command" )
					
add_custom_command ( TARGET ${GUI_NAME}
                    POST_BUILD
                    COMMAND echo TARGET_FILE:${GUI_NAME}=$<TARGET_FILE:${GUI_NAME}>
                    COMMAND echo TARGET_FILE_BASE_NAME:${GUI_NAME}=$<TARGET_FILE_BASE_NAME:${GUI_NAME}>
					COMMAND echo TARGET_FILE_PREFIX:${GUI_NAME}=$<TARGET_FILE_PREFIX:${GUI_NAME}>
					COMMAND echo TARGET_FILE_SUFFIX:${GUI_NAME}=$<TARGET_FILE_SUFFIX:${GUI_NAME}>
					COMMAND echo TARGET_FILE_NAME:${GUI_NAME}=$<TARGET_FILE_NAME:${GUI_NAME}>
                    COMMENT "Run custom command" )

include(InstallRequiredSystemLibraries)
message ( STATUS "CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS = ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}" )
message ( STATUS "CMAKE_INSTALL_RUNTIME_DESTINATION = ${CMAKE_INSTALL_RUNTIME_DESTINATION}" )

#install(
#    PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
#    DESTINATION ${CMAKE_INSTALL_RUNTIME_DESTINATION}
#)

install ( TARGETS ${GUI_NAME} ${PRJ_BIN_DIR_FILES}
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

set ( CPACK_PACKAGE_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_MAJOR}" )
set ( CPACK_PACKAGE_VERSION_MINOR "${${PROJECT_NAME}_VERSION_MINOR}" )
include(CPack)

install ( DIRECTORY "${OUTPUT_BIN_DIR}" DESTINATION "." )

