cmake_minimum_required(VERSION 3.30)

if(NOT DEFINED ONEFLOW_SOURCE_ROOT)
    message(FATAL_ERROR "ONEFLOW_SOURCE_ROOT is required")
endif()

set(_work "${CMAKE_CURRENT_BINARY_DIR}/source discovery fixture")
set(_source "${_work}/source")
set(_binary "${_work}/binary")
file(REMOVE_RECURSE "${_work}")
file(MAKE_DIRECTORY
    "${_source}/include"
    "${_source}/src/nested"
    "${_source}/legacy"
    "${_source}/test"
    "${_source}/CMakeFiles"
    "${_source}/build"
    "${_source}/empty"
)

file(WRITE "${_source}/include/Module.h" "#pragma once\n")
file(WRITE "${_source}/src/nested/NewStyle.cpp" "void new_style() {}\n")
file(WRITE "${_source}/legacy/Legacy.cpp" "void legacy_style() {}\n")
file(WRITE "${_source}/test/LegacyTestModule.cpp" "void legacy_test_module() {}\n")
file(WRITE "${_source}/build/ShouldNotCompile.cpp" "void generated_only() {}\n")
file(WRITE "${_source}/CMakeFiles/ShouldNotCompile.cpp" "void cmake_generated() {}\n")

set(_fixture_cmake "${_source}/CMakeLists.txt")
file(WRITE "${_fixture_cmake}" [=[
cmake_minimum_required(VERSION 3.30)
project(OneFLOWDiscoveryFixture NONE)
include("${ONEFLOW_DISCOVERY_MODULE}")
oneflow_discover_module_sources(
    "${CMAKE_CURRENT_SOURCE_DIR}" discovered_sources discovered_headers discovered_includes)
oneflow_discover_module_sources(
    "${CMAKE_CURRENT_SOURCE_DIR}/empty" empty_sources empty_headers empty_includes)
if(empty_sources OR empty_headers OR empty_includes)
    message(FATAL_ERROR "An empty module must produce empty discovery results")
endif()
file(WRITE "${CMAKE_BINARY_DIR}/discovery.txt"
    "SOURCES\n${discovered_sources}\nHEADERS\n${discovered_headers}\nINCLUDES\n${discovered_includes}\n")
]=])

function(configure_fixture)
    execute_process(
        COMMAND "${CMAKE_COMMAND}"
            -S "${_source}"
            -B "${_binary}"
            "-DONEFLOW_DISCOVERY_MODULE=${ONEFLOW_SOURCE_ROOT}/cmake/modules/OneFLOWSourceDiscovery.cmake"
        RESULT_VARIABLE _status
        OUTPUT_VARIABLE _stdout
        ERROR_VARIABLE _stderr
    )
    if(NOT _status EQUAL 0)
        message(FATAL_ERROR
            "Fixture configure failed (${_status}).\n${_stdout}\n${_stderr}")
    endif()
    file(READ "${_binary}/discovery.txt" _manifest)
    set(_manifest "${_manifest}" PARENT_SCOPE)
endfunction()

configure_fixture()
if(NOT _manifest MATCHES "Legacy[.]cpp" OR NOT _manifest MATCHES "NewStyle[.]cpp")
    message(FATAL_ERROR "Legacy and include/src source layouts were not both discovered:\n${_manifest}")
endif()
if(NOT _manifest MATCHES "LegacyTestModule[.]cpp")
    message(FATAL_ERROR "A valid legacy module under a directory named test was excluded:\n${_manifest}")
endif()
if(_manifest MATCHES "ShouldNotCompile[.]cpp")
    message(FATAL_ERROR "Generated build/CMakeFiles sources leaked into discovery:\n${_manifest}")
endif()
if(NOT _manifest MATCHES "Module[.]h")
    message(FATAL_ERROR "Module header was not discovered:\n${_manifest}")
endif()
if(NOT _manifest MATCHES "include")
    message(FATAL_ERROR "Header include directory was not discovered:\n${_manifest}")
endif()

# Verify that CONFIGURE_DEPENDS-compatible discovery sees additions and
# removals on a subsequent configure without touching the source repository.
file(WRITE "${_source}/src/nested/AddedLater.cpp" "void added_later() {}\n")
configure_fixture()
if(NOT _manifest MATCHES "AddedLater[.]cpp")
    message(FATAL_ERROR "A newly added source file was not discovered:\n${_manifest}")
endif()

file(REMOVE "${_source}/src/nested/AddedLater.cpp")
configure_fixture()
if(_manifest MATCHES "AddedLater[.]cpp")
    message(FATAL_ERROR "A removed source file remained in discovery:\n${_manifest}")
endif()

message(STATUS "OneFLOW CMake source discovery tests passed")
