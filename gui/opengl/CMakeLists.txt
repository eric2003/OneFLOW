project ( OneFLOW_OpenGL LANGUAGES C CXX )

set ( PRJ_SRC_LIST )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

#set ( MY_GLFW "d:/software/GLFW" )
#set ( MY_GLAD "d:/software/glad" )
set ( MY_GLFW "/home/kylinlu/work_gitlab/oneflow/glfw" )
set ( MY_GLAD "/home/kylinlu/work_gitlab/oneflow/glad" )

set ( MY_GLFW_INC "${MY_GLFW}/include" )
set ( MY_GLFW_LIB "${MY_GLFW}/build/src/libglfw.so" )
#set ( MY_GLFW_LIB "${MY_GLFW}/lib/glfw3dll.lib" )

set ( MY_GLAD_INC "${MY_GLAD}/include" )
set ( MY_GLAD_SRC "${MY_GLAD}/src" )

set ( glad_c_files "${MY_GLAD_SRC}/glad.c" )

message ( STATUS "CMAKE_SOURCE_DIR = ${CMAKE_SOURCE_DIR}" )
message ( STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}" )
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}" )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}" )


file ( GLOB root_header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h" )
file ( GLOB root_src_files "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" )
list ( APPEND root_src_files ${glad_c_files} )

list ( APPEND PRJ_INCLUDE_DIRS . )
list ( APPEND PRJ_INCLUDE_DIRS ${MY_GLAD_INC} )
list ( APPEND PRJ_INCLUDE_DIRS ${MY_GLFW_INC} )

list ( APPEND PRJ_SRC_LIST ${root_src_files} )
list ( APPEND PRJ_LIBRARIES ${MY_GLFW_LIB} )

add_executable ( ${PROJECT_NAME} ${PRJ_SRC_LIST} ${root_header_files} )

target_include_directories ( ${PROJECT_NAME}
	PRIVATE 
        ${PRJ_INCLUDE_DIRS}
)

target_link_libraries ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)
