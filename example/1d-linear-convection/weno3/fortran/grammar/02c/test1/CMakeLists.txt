
add_executable( main1
    main1.f90
	../src/sub1.f90
	../src/sub2.f90
)

target_compile_features ( main1
	PRIVATE 
		${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( main1
	PRIVATE
	   ${PRJ_COMPILE_DEFINITIONS} 
)