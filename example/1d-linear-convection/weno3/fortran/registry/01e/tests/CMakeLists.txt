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

message(STATUS "Tests configured")