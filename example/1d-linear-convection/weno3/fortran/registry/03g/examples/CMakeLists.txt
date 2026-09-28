# examples/CMakeLists.txt (修正版)
message(STATUS "配置示例程序...")

# 主示例程序：ENO/WENO对比（当前版本）
add_executable(run_eno_weno
    run_eno_weno.f90
)

target_link_libraries(run_eno_weno
    PRIVATE
        solver
        infrastructure
        core
        physics
        manager
        results
        boundary
        initial_condition
)

# 集成版本示例
add_executable(run_eno_weno_integrated
    run_eno_weno_integrated.f90
)

target_link_libraries(run_eno_weno_integrated
    PRIVATE
        solver_integrated  # 主要模块
        infrastructure
        base
        boundary          # 集成求解器需要的边界条件
        initial_condition # 集成求解器需要的初始条件
)

install(TARGETS run_eno_weno run_eno_weno_integrated
    RUNTIME DESTINATION bin/examples
)

message(STATUS "示例程序配置完成")