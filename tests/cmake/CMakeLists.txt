add_test(
    NAME CMakeSourceDiscovery
    COMMAND "\${CMAKE_COMMAND}"
        "-DONEFLOW_SOURCE_ROOT=\${PROJECT_SOURCE_DIR}"
        -P "\${CMAKE_CURRENT_SOURCE_DIR}/test_source_discovery.cmake"
)
set_tests_properties(CMakeSourceDiscovery PROPERTIES
    LABELS "cmake;infrastructure;unit"
)
