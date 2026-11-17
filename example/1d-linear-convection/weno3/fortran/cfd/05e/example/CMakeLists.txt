set(SOURCES
	run_eno_weno.f90
)

add_executable(${PROJECT_NAME} ${SOURCES})

if(MSVC)
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
endif()


target_link_libraries(${PROJECT_NAME}
    PRIVATE 
        cfd_lib
)

target_include_directories(${PROJECT_NAME} 
	PRIVATE 
		${CMAKE_Fortran_MODULE_DIRECTORY}
)
