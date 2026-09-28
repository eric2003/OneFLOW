# tests/CMakeLists.txt
message(STATUS "Configuring tests...")

add_executable(test_minimal_simple test_minimal_simple.f90)

message(STATUS "CMAKE_Fortran_MODULE_DIRECTORY=${CMAKE_Fortran_MODULE_DIRECTORY}")

#target_include_directories( test_minimal_simple
#	PRIVATE 
#		${CMAKE_Fortran_MODULE_DIRECTORY}
#)

target_link_libraries( test_minimal_simple
	PRIVATE 
		infrastructure
)

add_executable(test_simple_link test_simple_link.f90)
target_link_libraries(test_simple_link
    PRIVATE 
        reconstructor
        flux
)


#add_executable(test_factory_simple test_factory_simple.f90)

#target_link_libraries( test_factory_simple
#	PRIVATE 
#		core
#		infrastructure
#		reconstructor
#		flux
#)

add_executable(test_factory_simple test_factory_simple.f90)
target_link_libraries(test_factory_simple
    PRIVATE 
        core
        infrastructure
        reconstructor
        flux
)