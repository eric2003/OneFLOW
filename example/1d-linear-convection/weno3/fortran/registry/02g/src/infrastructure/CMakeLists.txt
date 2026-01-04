# src/infrastructure/CMakeLists.txt
message(STATUS "Configuring infrastructure module...")

add_library(infrastructure STATIC
    config.f90
    mesh.f90
    domain.f90      # 新增
    solution.f90    # 新增
)

target_link_libraries(infrastructure PRIVATE base)

set_target_properties(infrastructure PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "Infrastructure module configured")