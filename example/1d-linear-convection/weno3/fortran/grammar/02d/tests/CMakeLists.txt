
message(STATUS "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_BINARY_DIR/src/Debug/mylib.lib=${CMAKE_BINARY_DIR}/src/Debug/mylib.lib")


add_executable( ${PROJECT_NAME}
    main.f90
)

target_compile_features ( ${PROJECT_NAME} 
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)

target_link_libraries( ${PROJECT_NAME}
	PRIVATE 
		${CMAKE_BINARY_DIR}/src/Debug/mylib.lib
)