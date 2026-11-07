include_guard(GLOBAL)

# Discover implementation and header files below one module root.
#
# The first four arguments are output-variable names, matching the existing
# CMake helper style used by OneFLOW:
#   oneflow_discover_module_sources(root sources headers include_dirs)
#
# Both legacy layouts (implementation files directly below a module or
# arbitrary legacy subdirectories) and the include/src convention are
# accepted. Test, build and VCS directories are never production inputs.
function(oneflow_discover_module_sources module_root source_var header_var include_var)
    if(NOT IS_DIRECTORY "${module_root}")
        message(FATAL_ERROR
            "OneFLOW source discovery root does not exist: ${module_root}")
    endif()

    get_filename_component(_module_root "${module_root}" ABSOLUTE)
    file(GLOB_RECURSE _candidate_files
        CONFIGURE_DEPENDS
        LIST_DIRECTORIES FALSE
        "${_module_root}/*.cpp"
        "${_module_root}/*.h"
        "${_module_root}/*.hpp"
    )

    set(_sources)
    set(_headers)
    set(_include_dirs)

    foreach(_file IN LISTS _candidate_files)
        file(RELATIVE_PATH _relative "${_module_root}" "${_file}")
        file(TO_CMAKE_PATH "${_relative}" _relative)

        # Do not let test fixtures, generated build trees, or repository
        # metadata leak into production module discovery.
        if(_relative MATCHES "(^|/)(build|CMakeFiles|\\.git)(/|$)")
            continue()
        endif()

        get_filename_component(_parent "${_file}" DIRECTORY)
        get_filename_component(_parent_name "${_parent}" NAME)

        if(_file MATCHES "\\.(cpp)$")
            # Headers and implementation files have distinct ownership.
            if(_parent_name STREQUAL "include")
                continue()
            endif()
            list(APPEND _sources "${_file}")
        else()
            list(APPEND _headers "${_file}")
            list(APPEND _include_dirs "${_parent}")
        endif()
    endforeach()

    list(REMOVE_DUPLICATES _sources)
    list(REMOVE_DUPLICATES _headers)
    list(REMOVE_DUPLICATES _include_dirs)
    list(SORT _sources)
    list(SORT _headers)
    list(SORT _include_dirs)

    set(${source_var} "${_sources}" PARENT_SCOPE)
    set(${header_var} "${_headers}" PARENT_SCOPE)
    set(${include_var} "${_include_dirs}" PARENT_SCOPE)
endfunction()
