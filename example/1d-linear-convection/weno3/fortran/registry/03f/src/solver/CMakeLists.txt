# src/solver/CMakeLists.txt（统一管理所有求解器）
message(STATUS "配置求解器模块...")

# 基础求解器库
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
        results
        boundary        # 添加依赖
        initial_condition # 添加依赖
)

set_target_properties(solver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# 集成求解器库（在同一个目录下）
add_library(solver_integrated STATIC
    solver_integrated.f90
)

target_link_libraries(solver_integrated
    PRIVATE
        solver           # 基础求解器
        infrastructure
        base
        boundary
        initial_condition
        results          # 如果需要保存结果
)

set_target_properties(solver_integrated PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "求解器模块配置完成")