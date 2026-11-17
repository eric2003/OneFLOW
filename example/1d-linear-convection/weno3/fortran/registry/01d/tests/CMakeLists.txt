# tests/CMakeLists.txt
message(STATUS "Configuring tests...")

# Simple functionality test
add_executable(test_simple
    test_minimal_simple.f90
)
target_link_libraries(test_simple
    core
    infrastructure
)
target_include_directories(test_simple PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Factory pattern test
add_executable(test_factory
    test_factory.f90
)
target_link_libraries(test_factory
    core
    infrastructure
    reconstructor
    flux
)
target_include_directories(test_factory PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Test targets
add_custom_target(run_all_tests
    COMMAND echo "=== Running All Tests ==="
    COMMAND echo ""
    COMMAND echo "1. Simple functionality test..."
    COMMAND test_simple
    COMMAND echo ""
    COMMAND echo "2. Factory pattern test..."
    COMMAND test_factory
    COMMAND echo ""
    COMMAND echo "=== All Tests Completed ==="
    DEPENDS test_simple test_factory
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

message(STATUS "Tests configured")