cmake_minimum_required(VERSION 3.30)

if(NOT DEFINED ONEFLOW_SOURCE_ROOT)
    message(FATAL_ERROR "ONEFLOW_SOURCE_ROOT is required")
endif()

set(_work "${CMAKE_CURRENT_BINARY_DIR}/legacy discovery parity")
set(_source "${_work}/source")
set(_binary "${_work}/binary")
file(REMOVE_RECURSE "${_work}")
file(MAKE_DIRECTORY "${_source}")

file(WRITE "${_source}/CMakeLists.txt" [=[
cmake_minimum_required(VERSION 3.30)
project(OneFLOWLegacyDiscoveryParity NONE)

include("${ONEFLOW_SOURCE_ROOT}/cmake/UtilFunction.cmake")
include("${ONEFLOW_SOURCE_ROOT}/cmake/modules/OneFLOWSourceDiscovery.cmake")

ConstructSolutionDirTree(
    "${ONEFLOW_SOURCE_ROOT}/codes"
    legacy_headers legacy_sources legacy_includes)
oneflow_discover_module_sources(
    "${ONEFLOW_SOURCE_ROOT}/codes"
    discovered_sources discovered_headers discovered_includes)

list(SORT legacy_sources)
list(SORT legacy_headers)
list(SORT legacy_includes)

if(NOT legacy_sources STREQUAL discovered_sources)
    message(FATAL_ERROR
        "Source discovery differs from the legacy implementation.\n"
        "Legacy: ${legacy_sources}\nNew: ${discovered_sources}")
endif()
if(NOT legacy_headers STREQUAL discovered_headers)
    message(FATAL_ERROR
        "Header discovery differs from the legacy implementation.\n"
        "Legacy: ${legacy_headers}\nNew: ${discovered_headers}")
endif()
if(NOT legacy_includes STREQUAL discovered_includes)
    message(FATAL_ERROR
        "Include-directory discovery differs from the legacy implementation.\n"
        "Legacy: ${legacy_includes}\nNew: ${discovered_includes}")
endif()

message(STATUS "Legacy and new OneFLOW source discovery results match")
]=])

execute_process(
    COMMAND "${CMAKE_COMMAND}"
        -S "${_source}"
        -B "${_binary}"
        "-DONEFLOW_SOURCE_ROOT=${ONEFLOW_SOURCE_ROOT}"
    RESULT_VARIABLE _status
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr
)
if(NOT _status EQUAL 0)
    message(FATAL_ERROR
        "Legacy discovery parity check failed (${_status}).\n${_stdout}\n${_stderr}")
endif()

message(STATUS "OneFLOW legacy discovery parity test passed")
