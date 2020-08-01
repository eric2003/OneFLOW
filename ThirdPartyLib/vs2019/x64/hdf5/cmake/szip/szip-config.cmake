#-----------------------------------------------------------------------------
# SZIP Config file for compiling against SZIP build directory
#-----------------------------------------------------------------------------

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was szip-config.cmake.in                            ########

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

string(TOUPPER szip SZIP_PACKAGE_NAME)

set (${SZIP_PACKAGE_NAME}_VALID_COMPONENTS static shared)

#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
set (${SZIP_PACKAGE_NAME}_ENABLE_ENCODING      ON)
set (${SZIP_PACKAGE_NAME}_BUILD_SHARED_LIBS    ON)
set (${SZIP_PACKAGE_NAME}_EXPORT_LIBRARIES     szip-static;szip-shared)

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
set (${SZIP_PACKAGE_NAME}_INCLUDE_DIR "${PACKAGE_PREFIX_DIR}/include")

set (${SZIP_PACKAGE_NAME}_SHARE_DIR "${PACKAGE_PREFIX_DIR}/cmake")
set_and_check (${SZIP_PACKAGE_NAME}_BUILD_DIR "${PACKAGE_PREFIX_DIR}")

#-----------------------------------------------------------------------------
# Version Strings
#-----------------------------------------------------------------------------
set (${SZIP_PACKAGE_NAME}_VERSION_STRING 2.1)
set (${SZIP_PACKAGE_NAME}_VERSION_MAJOR  2.1)
set (${SZIP_PACKAGE_NAME}_VERSION_MINOR  1)

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already build SZIP as a subproject
#-----------------------------------------------------------------------------
if (NOT TARGET "szip")
  include (${PACKAGE_PREFIX_DIR}/cmake/szip/szip-targets.cmake)
endif ()

# Handle default component(static) :
if (NOT ${SZIP_PACKAGE_NAME}_FIND_COMPONENTS)
    set (${SZIP_PACKAGE_NAME}_FIND_COMPONENTS static)
    set (${SZIP_PACKAGE_NAME}_FIND_REQUIRED_static true)
endif ()

# Handle requested components:
list (REMOVE_DUPLICATES ${SZIP_PACKAGE_NAME}_FIND_COMPONENTS)
foreach (comp IN LISTS ${SZIP_PACKAGE_NAME}_FIND_COMPONENTS)
    list (FIND ${SZIP_PACKAGE_NAME}_EXPORT_LIBRARIES "szip-${comp}" HAVE_COMP) 
    if (${HAVE_COMP} LESS 0) 
      set (${SZIP_PACKAGE_NAME}_${comp}_FOUND 0)
    else ()
      set (${SZIP_PACKAGE_NAME}_${comp}_FOUND 1)
      string(TOUPPER ${SZIP_PACKAGE_NAME}_${comp}_LIBRARY COMP_LIBRARY)
      set (${COMP_LIBRARY} ${${COMP_LIBRARY}} szip-${comp})
    endif ()
endforeach ()

check_required_components (${SZIP_PACKAGE_NAME})
