cmake_minimum_required ( VERSION 3.26 )

project ( testprj )

set ( PRJ_INCLUDE_DIRS )
set ( PRJ_COMPILE_FEATURES )
set ( PRJ_LIBRARIES )

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

find_package(Eigen3 CONFIG REQUIRED)

message ( STATUS "EIGEN3_INCLUDE_DIR    = ${EIGEN3_INCLUDE_DIR}" )

list ( APPEND PRJ_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )
list ( APPEND PRJ_LIBRARIES Eigen3::Eigen )

set_property ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )

add_executable ( ${PROJECT_NAME}
    main.cpp
)

target_include_directories ( ${PROJECT_NAME}
    PRIVATE 
        ${PRJ_INCLUDE_DIRS}
)

target_link_libraries( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
    PRIVATE 
        ${PRJ_COMPILE_FEATURES}
)