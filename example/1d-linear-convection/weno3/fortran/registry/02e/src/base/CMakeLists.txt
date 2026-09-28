# src/base/CMakeLists.txt
message(STATUS "Configuring base module...")

add_library(base STATIC
    modules.f90
)

set_target_properties(base PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "Base module configured")

# src/core/CMakeLists.txt  
message(STATUS "Configuring core module...")