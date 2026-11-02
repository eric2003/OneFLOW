cmake_minimum_required ( VERSION 3.15 )

project ( vtk_prj )

set ( PRJ_SRC_LIST )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

find_package ( VTK )

if ( VTK_FOUND )
	message ( STATUS " VTK_FOUND =  ${VTK_FOUND}" )
	message ( STATUS " VTK_MAJOR_VERSION = ${VTK_MAJOR_VERSION}" )
	message ( STATUS " VTK_MINOR_VERSION = ${VTK_MINOR_VERSION}" )
	message ( STATUS " VTK_BUILD_VERSION = ${VTK_BUILD_VERSION}" )
	message ( STATUS " VTK_LIBRARY_DIRS  = ${VTK_LIBRARY_DIRS}" )
	message ( STATUS " VTK_LIBRARIES  = ${VTK_LIBRARIES}" )
endif()

file ( GLOB root_header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" )
file ( GLOB root_src_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" )

list ( APPEND PRJ_INCLUDE_DIRS . )
list ( APPEND PRJ_INCLUDE_DIRS ${VTK_USE_FILE} )

list ( APPEND PRJ_SRC_LIST ${root_src_files} )
list ( APPEND PRJ_LIBRARIES ${VTK_LIBRARIES} )

add_executable ( ${PROJECT_NAME} ${PRJ_SRC_LIST} )

target_include_directories ( ${PROJECT_NAME}
	PRIVATE 
        ${PRJ_INCLUDE_DIRS}
)

target_link_libraries ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)

# vtk_module_autoinit is needed
vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
)
