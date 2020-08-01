#-----------------------------------------------------------------------------
# ZLIB Config file for compiling against ZLIB build directory
#-----------------------------------------------------------------------------

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was zlib-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

string(TOUPPER zlib ZLIB_PACKAGE_NAME)

set (${ZLIB_PACKAGE_NAME}_VALID_COMPONENTS static shared)

#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
set (${ZLIB_PACKAGE_NAME}_BUILD_SHARED_LIBS    ON)
set (${ZLIB_PACKAGE_NAME}_EXPORT_LIBRARIES     zlib-static;zlib-shared)

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
set (${ZLIB_PACKAGE_NAME}_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")

set (${ZLIB_PACKAGE_NAME}_SHARE_DIR "${PACKAGE_PREFIX_DIR}/cmake")
set_and_check (${ZLIB_PACKAGE_NAME}_BUILD_DIR "${PACKAGE_PREFIX_DIR}")

#-----------------------------------------------------------------------------
# Version Strings
#-----------------------------------------------------------------------------
set (${ZLIB_PACKAGE_NAME}_VERSION_STRING 1.2)
set (${ZLIB_PACKAGE_NAME}_VERSION_MAJOR  1.2)
set (${ZLIB_PACKAGE_NAME}_VERSION_MINOR  11)

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already build ZLIB as a subproject
#-----------------------------------------------------------------------------
if (NOT TARGET "zlib")
  include (${PACKAGE_PREFIX_DIR}/cmake/zlib/zlib-targets.cmake)
endif ()

# Handle default component(static) :
if (NOT ${ZLIB_PACKAGE_NAME}_FIND_COMPONENTS)
    set (${ZLIB_PACKAGE_NAME}_FIND_COMPONENTS static)
    set (${ZLIB_PACKAGE_NAME}_FIND_REQUIRED_static true)
endif ()

# Handle requested components:
list (REMOVE_DUPLICATES ${ZLIB_PACKAGE_NAME}_FIND_COMPONENTS)
foreach (comp IN LISTS ${ZLIB_PACKAGE_NAME}_FIND_COMPONENTS)
    list (FIND ${ZLIB_PACKAGE_NAME}_EXPORT_LIBRARIES "zlib-${comp}" HAVE_COMP) 
    if (${HAVE_COMP} LESS 0) 
      set (${ZLIB_PACKAGE_NAME}_${comp}_FOUND 0)
    else ()
      set (${ZLIB_PACKAGE_NAME}_${comp}_FOUND 1)
      string(TOUPPER ${ZLIB_PACKAGE_NAME}_${comp}_LIBRARY COMP_LIBRARY)
      set (${COMP_LIBRARY} ${${COMP_LIBRARY}} zlib-${comp})
    endif ()
endforeach ()

check_required_components (${ZLIB_PACKAGE_NAME})
