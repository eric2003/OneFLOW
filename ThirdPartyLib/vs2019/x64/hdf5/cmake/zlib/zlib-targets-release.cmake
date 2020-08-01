#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hdf5::zlib-static" for configuration "Release"
set_property(TARGET hdf5::zlib-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::zlib-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libzlib.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::zlib-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::zlib-static "${_IMPORT_PREFIX}/lib/libzlib.lib" )

# Import target "hdf5::zlib-shared" for configuration "Release"
set_property(TARGET hdf5::zlib-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::zlib-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/zlib.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/zlib.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::zlib-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::zlib-shared "${_IMPORT_PREFIX}/lib/zlib.lib" "${_IMPORT_PREFIX}/bin/zlib.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
