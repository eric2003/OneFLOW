#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hdf5::hdf5-static" for configuration "Release"
set_property(TARGET hdf5::hdf5-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5-static "${_IMPORT_PREFIX}/lib/libhdf5.lib" )

# Import target "hdf5::hdf5-shared" for configuration "Release"
set_property(TARGET hdf5::hdf5-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/hdf5.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/hdf5.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5-shared "${_IMPORT_PREFIX}/lib/hdf5.lib" "${_IMPORT_PREFIX}/bin/hdf5.dll" )

# Import target "hdf5::hdf5_tools-static" for configuration "Release"
set_property(TARGET hdf5::hdf5_tools-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_tools-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_tools.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_tools-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_tools-static "${_IMPORT_PREFIX}/lib/libhdf5_tools.lib" )

# Import target "hdf5::hdf5_tools-shared" for configuration "Release"
set_property(TARGET hdf5::hdf5_tools-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_tools-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/hdf5_tools.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/hdf5_tools.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_tools-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_tools-shared "${_IMPORT_PREFIX}/lib/hdf5_tools.lib" "${_IMPORT_PREFIX}/bin/hdf5_tools.dll" )

# Import target "hdf5::h5diff" for configuration "Release"
set_property(TARGET hdf5::h5diff APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5diff PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5diff.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5diff )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5diff "${_IMPORT_PREFIX}/bin/h5diff.exe" )

# Import target "hdf5::h5diff-shared" for configuration "Release"
set_property(TARGET hdf5::h5diff-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5diff-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5diff-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5diff-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5diff-shared "${_IMPORT_PREFIX}/bin/h5diff-shared.exe" )

# Import target "hdf5::h5ls" for configuration "Release"
set_property(TARGET hdf5::h5ls APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5ls PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5ls.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5ls )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5ls "${_IMPORT_PREFIX}/bin/h5ls.exe" )

# Import target "hdf5::h5ls-shared" for configuration "Release"
set_property(TARGET hdf5::h5ls-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5ls-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5ls-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5ls-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5ls-shared "${_IMPORT_PREFIX}/bin/h5ls-shared.exe" )

# Import target "hdf5::h5debug" for configuration "Release"
set_property(TARGET hdf5::h5debug APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5debug PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5debug.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5debug )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5debug "${_IMPORT_PREFIX}/bin/h5debug.exe" )

# Import target "hdf5::h5repart" for configuration "Release"
set_property(TARGET hdf5::h5repart APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5repart PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5repart.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5repart )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5repart "${_IMPORT_PREFIX}/bin/h5repart.exe" )

# Import target "hdf5::h5mkgrp" for configuration "Release"
set_property(TARGET hdf5::h5mkgrp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5mkgrp PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5mkgrp.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5mkgrp )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5mkgrp "${_IMPORT_PREFIX}/bin/h5mkgrp.exe" )

# Import target "hdf5::h5import" for configuration "Release"
set_property(TARGET hdf5::h5import APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5import PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5import.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5import )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5import "${_IMPORT_PREFIX}/bin/h5import.exe" )

# Import target "hdf5::h5repack" for configuration "Release"
set_property(TARGET hdf5::h5repack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5repack PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5repack.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5repack )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5repack "${_IMPORT_PREFIX}/bin/h5repack.exe" )

# Import target "hdf5::h5repack-shared" for configuration "Release"
set_property(TARGET hdf5::h5repack-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5repack-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5repack-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5repack-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5repack-shared "${_IMPORT_PREFIX}/bin/h5repack-shared.exe" )

# Import target "hdf5::h5jam" for configuration "Release"
set_property(TARGET hdf5::h5jam APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5jam PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5jam.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5jam )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5jam "${_IMPORT_PREFIX}/bin/h5jam.exe" )

# Import target "hdf5::h5unjam" for configuration "Release"
set_property(TARGET hdf5::h5unjam APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5unjam PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5unjam.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5unjam )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5unjam "${_IMPORT_PREFIX}/bin/h5unjam.exe" )

# Import target "hdf5::h5copy" for configuration "Release"
set_property(TARGET hdf5::h5copy APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5copy PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5copy.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5copy )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5copy "${_IMPORT_PREFIX}/bin/h5copy.exe" )

# Import target "hdf5::h5copy-shared" for configuration "Release"
set_property(TARGET hdf5::h5copy-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5copy-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5copy-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5copy-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5copy-shared "${_IMPORT_PREFIX}/bin/h5copy-shared.exe" )

# Import target "hdf5::h5stat" for configuration "Release"
set_property(TARGET hdf5::h5stat APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5stat PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5stat.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5stat )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5stat "${_IMPORT_PREFIX}/bin/h5stat.exe" )

# Import target "hdf5::h5stat-shared" for configuration "Release"
set_property(TARGET hdf5::h5stat-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5stat-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5stat-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5stat-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5stat-shared "${_IMPORT_PREFIX}/bin/h5stat-shared.exe" )

# Import target "hdf5::h5dump" for configuration "Release"
set_property(TARGET hdf5::h5dump APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5dump PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5dump.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5dump )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5dump "${_IMPORT_PREFIX}/bin/h5dump.exe" )

# Import target "hdf5::h5dump-shared" for configuration "Release"
set_property(TARGET hdf5::h5dump-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h5dump-shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h5dump-shared.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h5dump-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h5dump-shared "${_IMPORT_PREFIX}/bin/h5dump-shared.exe" )

# Import target "hdf5::hdf5_hl-static" for configuration "Release"
set_property(TARGET hdf5::hdf5_hl-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_hl-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_hl.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_hl-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_hl-static "${_IMPORT_PREFIX}/lib/libhdf5_hl.lib" )

# Import target "hdf5::hdf5_hl-shared" for configuration "Release"
set_property(TARGET hdf5::hdf5_hl-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_hl-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/hdf5_hl.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/hdf5_hl.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_hl-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_hl-shared "${_IMPORT_PREFIX}/lib/hdf5_hl.lib" "${_IMPORT_PREFIX}/bin/hdf5_hl.dll" )

# Import target "hdf5::gif2h5" for configuration "Release"
set_property(TARGET hdf5::gif2h5 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::gif2h5 PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/gif2h5.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::gif2h5 )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::gif2h5 "${_IMPORT_PREFIX}/bin/gif2h5.exe" )

# Import target "hdf5::h52gif" for configuration "Release"
set_property(TARGET hdf5::h52gif APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::h52gif PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/h52gif.exe"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::h52gif )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::h52gif "${_IMPORT_PREFIX}/bin/h52gif.exe" )

# Import target "hdf5::hdf5_cpp-static" for configuration "Release"
set_property(TARGET hdf5::hdf5_cpp-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_cpp-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_cpp.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_cpp-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_cpp-static "${_IMPORT_PREFIX}/lib/libhdf5_cpp.lib" )

# Import target "hdf5::hdf5_cpp-shared" for configuration "Release"
set_property(TARGET hdf5::hdf5_cpp-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_cpp-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/hdf5_cpp.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/hdf5_cpp.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_cpp-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_cpp-shared "${_IMPORT_PREFIX}/lib/hdf5_cpp.lib" "${_IMPORT_PREFIX}/bin/hdf5_cpp.dll" )

# Import target "hdf5::hdf5_hl_cpp-static" for configuration "Release"
set_property(TARGET hdf5::hdf5_hl_cpp-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_hl_cpp-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_hl_cpp.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_hl_cpp-static )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_hl_cpp-static "${_IMPORT_PREFIX}/lib/libhdf5_hl_cpp.lib" )

# Import target "hdf5::hdf5_hl_cpp-shared" for configuration "Release"
set_property(TARGET hdf5::hdf5_hl_cpp-shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5::hdf5_hl_cpp-shared PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/hdf5_hl_cpp.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/hdf5_hl_cpp.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5::hdf5_hl_cpp-shared )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5::hdf5_hl_cpp-shared "${_IMPORT_PREFIX}/lib/hdf5_hl_cpp.lib" "${_IMPORT_PREFIX}/bin/hdf5_hl_cpp.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
