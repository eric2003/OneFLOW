function ( BuildSlnTree currdir my_head_list my_src_list my_include_dirs )
    set ( tmp_header_list  "" )
    set ( tmp_src_list     "" )
    set ( tmp_include_dirs "" )

    message ( STATUS "BuildSlnTree: The currdir is ${currdir}" )
    file ( GLOB_RECURSE FOUND_FILES LIST_DIRECTORIES true RELATIVE ${currdir} * )

    #message( STATUS "BuildSlnTree: FOUND_FILES = ${FOUND_FILES}" )

    foreach ( child ${FOUND_FILES} )
        #message ( STATUS "BuildSlnTree:child=${child}" )
        set( candidate_dir ${currdir}/${child} )
        if ( IS_DIRECTORY ${candidate_dir} )
            set( local_dir ${candidate_dir} )
            set( group_name "${child}" )
        else()
            set( local_dir ${currdir} )
            #set( group_name "ROOT" )
            set( group_name  )
        endif()
        
        #message ( STATUS "BuildSlnTree:local_dir=${local_dir}, group_name=${group_name}" )
        
        if ( TRUE )
            file ( GLOB header_files "${local_dir}/*.h" )
            file ( GLOB src_files    "${local_dir}/*.cpp" )
            file ( GLOB cuda_src_files "${local_dir}/*.cu" )
            file ( GLOB hpp_files    "${local_dir}/*.hpp" )
            
            list ( APPEND header_files ${hpp_files} )
            
            list ( LENGTH header_files n_header_files )
            if ( NOT ( ${n_header_files} EQUAL 0 ) )
                list ( APPEND tmp_include_dirs ${local_dir} )
            endif()
            
            source_group ( "${group_name}" FILES ${header_files} )
            source_group ( "${group_name}" FILES ${src_files}    )
            source_group ( "${group_name}" FILES ${cuda_src_files})
            
            list ( APPEND tmp_header_list ${header_files} )
            list ( APPEND tmp_src_list  ${src_files} )
            list ( APPEND tmp_src_list  ${cuda_src_files} )
        endif()
    endforeach()
    set ( ${my_head_list} ${tmp_header_list} PARENT_SCOPE )
    set ( ${my_src_list} ${tmp_src_list} PARENT_SCOPE )
    set ( ${my_include_dirs} ${tmp_include_dirs} PARENT_SCOPE ) 
endfunction()