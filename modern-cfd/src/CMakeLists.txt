set ( SRC_ROOT ${CMAKE_CURRENT_SOURCE_DIR} )
list ( APPEND PRJ_INCLUDE_DIRS ${SRC_ROOT} )

message ( STATUS "SRC_ROOT=${SRC_ROOT}" )

BuildSlnTree( ${SRC_ROOT} PRJ_HEADER_LIST PRJ_SRC_LIST PRJ_LOCAL_INCLUDE_DIRS )
list ( APPEND PRJ_INCLUDE_DIRS ${PRJ_LOCAL_INCLUDE_DIRS} )

add_executable ( ${PROJECT_NAME}
    ${PRJ_HEADER_LIST}
    ${PRJ_SRC_LIST}
)

set_target_properties( ${PROJECT_NAME}
    PROPERTIES 
        POSITION_INDEPENDENT_CODE ON
)

if ( CUDAToolkit_FOUND )
set_target_properties ( ${PROJECT_NAME}
    PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)
endif()

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

target_compile_definitions ( ${PROJECT_NAME}
    PRIVATE
       ${PRJ_COMPILE_DEFINITIONS} 
)

install ( TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
