# examples/CMakeLists.txt
message(STATUS "配置示例程序...")

# 主示例程序：ENO/WENO对比
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
        results  # ← 新增链接结果模块
)

install(TARGETS run_eno_weno
    RUNTIME DESTINATION bin/examples
)

message(STATUS "示例程序配置完成")