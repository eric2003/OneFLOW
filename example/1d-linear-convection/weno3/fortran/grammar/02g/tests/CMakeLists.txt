
add_executable( ${PROJECT_NAME}
    main.f90
)

target_include_directories( ${PROJECT_NAME} 
	PRIVATE 
		${CMAKE_BINARY_DIR}/modules
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
		mylib
)