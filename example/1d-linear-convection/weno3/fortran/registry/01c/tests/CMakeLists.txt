# tests/CMakeLists.txt
message(STATUS "Configuring tests...")

# Minimal functionality test
add_executable(test_minimal
    test_minimal.f90
)
target_link_libraries(test_minimal
    core
    infrastructure
)
target_include_directories(test_minimal PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Test target
add_custom_target(run_tests
    COMMAND echo "=== Running Tests ==="
    COMMAND echo "Running minimal test..."
    COMMAND test_minimal
    COMMAND echo ""
    COMMAND echo "=== Tests Completed ==="
    DEPENDS test_minimal
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

message(STATUS "Tests configured")