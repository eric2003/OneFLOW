# examples/CMakeLists.txt
message(STATUS "配置示例程序...")

# 主示例程序：ENO/WENO对比
add_executable(example_eno_weno_comparison
    run_eno_weno.f90
)

target_link_libraries(example_eno_weno_comparison
    PRIVATE 
        solver
        infrastructure
        core
        physics
        manager
)

install(TARGETS example_eno_weno_comparison
    RUNTIME DESTINATION bin/examples
)

message(STATUS "示例程序配置完成")