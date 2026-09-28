# src/solver/CMakeLists.txt
message(STATUS "配置求解器模块...")

add_library(solver STATIC
    base.f90
    physics_solver.f90
)

target_link_libraries(solver
    PRIVATE
        infrastructure
        core
        physics
        manager
        results  # ← 新增链接结果模块
)

set_target_properties(solver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "求解器模块配置完成")