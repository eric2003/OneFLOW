function ( ConstructSolutionDirTree currdir my_head_list my_src_list my_include_dirs )
    set ( tmp_header_list  "" )
	set ( tmp_src_list     "" )
	set ( tmp_include_dirs "" )

	#message ( STATUS "The currdir is ${currdir}" )
	file ( GLOB_RECURSE FOUND_FILES LIST_DIRECTORIES true RELATIVE ${currdir} * )

	#message( STATUS "Files are ${FOUND_FILES}" )

    foreach ( child ${FOUND_FILES} )
		set( candidate_dir ${currdir}/${child} )
        if ( IS_DIRECTORY ${candidate_dir} )
			#message ( STATUS "The ${candidate_dir} is DIRECTORY" )
			file ( GLOB header_files "${candidate_dir}/*.h" )
			file ( GLOB src_files    "${candidate_dir}/*.cpp" )
			#file ( GLOB cuda_src_files "${candidate_dir}/*.cu" )
			file ( GLOB hpp_files    "${candidate_dir}/*.hpp" )
			
			list ( APPEND header_files ${hpp_files} )
			
			list ( LENGTH header_files n_header_files )
			#message ( STATUS "The n_header_files is ${n_header_files}" )
			if ( NOT ( ${n_header_files} EQUAL 0 ) )
				list ( APPEND tmp_include_dirs ${candidate_dir} )
			endif()
			
			source_group ( "${child}" FILES ${header_files} )
			source_group ( "${child}" FILES ${src_files}    )
			#source_group ( "${child}" FILES ${cuda_src_files}    )
			
			list ( APPEND tmp_header_list ${header_files} )
			list ( APPEND tmp_src_list  ${src_files} )
			#list ( APPEND tmp_src_list  ${cuda_src_files} )
		
			#message ( STATUS "The header_files is ${header_files}" )
			#message ( STATUS "The src_files is ${src_files}" )
			#message ( STATUS "The cuda_src_files is ${cuda_src_files}" )
        endif()
    endforeach()
    set ( ${my_head_list} ${tmp_header_list} PARENT_SCOPE )
    set ( ${my_src_list} ${tmp_src_list} PARENT_SCOPE )
	set ( ${my_include_dirs} ${tmp_include_dirs} PARENT_SCOPE )	
endfunction()

function ( GetCUDAFiles currdir my_head_list my_src_list my_include_dirs )
    set ( tmp_header_list  "" )
	set ( tmp_src_list     "" )
	set ( tmp_include_dirs "" )

	#message ( STATUS "The currdir is ${currdir}" )
	file ( GLOB_RECURSE FOUND_FILES LIST_DIRECTORIES true RELATIVE ${currdir} * )

    foreach ( child ${FOUND_FILES} )
		set( candidate_dir ${currdir}/${child} )
        if ( IS_DIRECTORY ${candidate_dir} )
			#message ( STATUS "The ${candidate_dir} is DIRECTORY" )
			file ( GLOB header_files "${candidate_dir}/*.h" )
			file ( GLOB cuda_src_files "${candidate_dir}/*.cu" )
				
			list ( LENGTH header_files n_header_files )
			if ( NOT ( ${n_header_files} EQUAL 0 ) )
				list ( APPEND tmp_header_list ${header_files} )
				list ( APPEND tmp_include_dirs ${candidate_dir} )				
			endif()
			list ( LENGTH cuda_src_files n_cuda_files )			
			if ( NOT ( ${n_cuda_files} EQUAL 0 ) )
				#source_group ( "${child}" FILES ${header_files} )
				source_group ( "${child}" FILES ${cuda_src_files}    )
				list ( APPEND tmp_src_list  ${cuda_src_files} )
				#message ( STATUS "The header_files is ${header_files}" )
				#message ( STATUS "The cuda_src_files is ${cuda_src_files}" )
			endif()
			
        endif()
    endforeach()
    set ( ${my_head_list} ${tmp_header_list} PARENT_SCOPE )
    set ( ${my_src_list} ${tmp_src_list} PARENT_SCOPE )
	set ( ${my_include_dirs} ${tmp_include_dirs} PARENT_SCOPE )	
endfunction()

function ( GetHIPFiles currdir my_head_list my_src_list my_include_dirs )
    set ( tmp_header_list  "" )
    set ( tmp_src_list     "" )
    set ( tmp_include_dirs "" )

    file ( GLOB_RECURSE FOUND_FILES LIST_DIRECTORIES true RELATIVE ${currdir} * )

    foreach ( child ${FOUND_FILES} )
        set( candidate_dir ${currdir}/${child} )
        if ( IS_DIRECTORY ${candidate_dir} )
            file ( GLOB header_files "${candidate_dir}/*.h" )
            file ( GLOB hip_src_files "${candidate_dir}/*.hip" )

            list ( LENGTH header_files n_header_files )
            if ( NOT ( ${n_header_files} EQUAL 0 ) )
                list ( APPEND tmp_header_list ${header_files} )
                list ( APPEND tmp_include_dirs ${candidate_dir} )
            endif()

            list ( LENGTH hip_src_files n_hip_files )
            if ( NOT ( ${n_hip_files} EQUAL 0 ) )
                source_group ( "${child}" FILES ${hip_src_files} )
                list ( APPEND tmp_src_list ${hip_src_files} )
            endif()
        endif()
    endforeach()

    set ( ${my_head_list} ${tmp_header_list} PARENT_SCOPE )
    set ( ${my_src_list} ${tmp_src_list} PARENT_SCOPE )
    set ( ${my_include_dirs} ${tmp_include_dirs} PARENT_SCOPE )
endfunction()


function ( DetectHIPArchitectures output_var )
    if ( CMAKE_HIP_ARCHITECTURES )
        set ( ${output_var} "${CMAKE_HIP_ARCHITECTURES}" PARENT_SCOPE )
        return()
    endif()

    if ( DEFINED ENV{ONEFLOW_HIP_ARCHITECTURES}
         AND NOT "$ENV{ONEFLOW_HIP_ARCHITECTURES}" STREQUAL "" )
        string ( REPLACE "," ";" detected_architectures
                 "$ENV{ONEFLOW_HIP_ARCHITECTURES}" )
        set ( ${output_var} "${detected_architectures}" PARENT_SCOPE )
        return()
    endif()

    # Prefer the architecture reported by an actually visible accelerator.
    # rocm_agent_enumerator lists compiler-supported targets and can contain
    # several unrelated gfx generations, so it is only a conservative fallback.
    find_program ( ONEFLOW_ROCMINFO rocminfo )
    if ( ONEFLOW_ROCMINFO )
        execute_process (
            COMMAND ${ONEFLOW_ROCMINFO}
            RESULT_VARIABLE rocminfo_status
            OUTPUT_VARIABLE rocminfo_output
            ERROR_QUIET
        )
        if ( rocminfo_status EQUAL 0 )
            string ( REGEX MATCHALL "gfx[0-9a-fA-F]+" detected_architectures
                     "${rocminfo_output}" )
            list ( REMOVE_ITEM detected_architectures gfx000 )
            list ( REMOVE_DUPLICATES detected_architectures )
            if ( detected_architectures )
                set ( ${output_var} "${detected_architectures}" PARENT_SCOPE )
                return()
            endif()
        endif()
    endif()

    find_program ( ONEFLOW_ROCM_AGENT_ENUMERATOR rocm_agent_enumerator )
    if ( ONEFLOW_ROCM_AGENT_ENUMERATOR )
        execute_process (
            COMMAND ${ONEFLOW_ROCM_AGENT_ENUMERATOR}
            RESULT_VARIABLE detector_status
            OUTPUT_VARIABLE detector_output
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if ( detector_status EQUAL 0 )
            string ( REGEX MATCHALL "gfx[0-9a-fA-F]+" detected_architectures
                     "${detector_output}" )
            list ( REMOVE_ITEM detected_architectures gfx000 )
            list ( REMOVE_DUPLICATES detected_architectures )
            list ( LENGTH detected_architectures detected_count )
            if ( detected_count EQUAL 1 )
                set ( ${output_var} "${detected_architectures}" PARENT_SCOPE )
                return()
            endif()
        endif()
    endif()

    set ( ${output_var} "" PARENT_SCOPE )
endfunction()

function ( AppendGlobalValue global_property value )
	#message ( STATUS "global_property = ${global_property}" )
	#message ( STATUS "value = ${value}" )
	get_property ( _local_var GLOBAL PROPERTY ${global_property} )
	
	list ( APPEND _local_var ${value} )
	set_property ( GLOBAL PROPERTY ${global_property} ${_local_var} )
endfunction()

function ( GetGlobalValue global_property value )
	#message ( STATUS "global_property = ${global_property}" )
	#message ( STATUS "value = ${value}" )
	get_property ( _local_var GLOBAL PROPERTY ${global_property} )
	set ( ${value} ${_local_var} PARENT_SCOPE )
endfunction()
