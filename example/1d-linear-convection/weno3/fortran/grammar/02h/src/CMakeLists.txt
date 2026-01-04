
set(mylibname "mylib" )

add_library(${mylibname} STATIC 
    sub1.f90
	sub2.f90
)

target_include_directories( ${mylibname}
	PRIVATE 
		${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_compile_features ( ${mylibname}
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${mylibname}
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)