# src/core/CMakeLists.txt  
message(STATUS "Configuring core module...")

add_library(core STATIC
    registry.f90
)

target_link_libraries(core PRIVATE base)

set_target_properties(core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "Core module configured")