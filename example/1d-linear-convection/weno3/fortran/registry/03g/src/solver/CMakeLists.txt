# src/solver/CMakeLists.txt（简化版）
message(STATUS "配置求解器模块...")

# 基础求解器库
add_library(solver STATIC
    base.f90
    physics_solver.f90
    residual_simple.f90  # 添加简化残差计算器
)

target_link_libraries(solver
    PRIVATE
        infrastructure
        core
        physics
        manager
        results
        boundary
        initial_condition
)

set_target_properties(solver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# 简单集成求解器库（演示用）
add_library(solver_integrated STATIC
    solver_integrated.f90
)

target_link_libraries(solver_integrated
    PRIVATE
        solver
        infrastructure
        base
        boundary
        initial_condition
)

set_target_properties(solver_integrated PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "求解器模块配置完成")