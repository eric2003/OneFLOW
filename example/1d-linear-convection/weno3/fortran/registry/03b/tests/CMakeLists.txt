# tests/CMakeLists.txt
message(STATUS "Configuring tests...")

# 基础设施测试
add_executable(test_infrastructure test_infrastructure.f90)
target_link_libraries(test_infrastructure
    PRIVATE 
        infrastructure
        core
)

# 注册系统测试
add_executable(test_registry test_registry.f90)
target_link_libraries(test_registry
    PRIVATE 
        core
        infrastructure
)

# 物理模块测试
add_executable(test_physics test_physics.f90)
target_link_libraries(test_physics
    PRIVATE 
        physics
        base
)

# 组件管理器测试
add_executable(test_component_manager test_component_manager.f90)
target_link_libraries(test_component_manager
    PRIVATE 
        manager
        infrastructure
)

# 配置物理测试
add_executable(test_config_physics test_config_physics.f90)
target_link_libraries(test_config_physics
    PRIVATE 
        infrastructure
        core
)

# 求解器基础测试
add_executable(test_solver_base test_solver_base.f90)
target_link_libraries(test_solver_base
    PRIVATE
        solver
        infrastructure
        core
)

# 物理求解器测试
add_executable(test_physics_solver test_physics_solver.f90)
target_link_libraries(test_physics_solver
    PRIVATE
        solver
        infrastructure
        core
        physics
        manager
)

add_executable(test_simple_link test_simple_link.f90)
target_link_libraries(test_simple_link
    PRIVATE 
        reconstructor
        flux
)


add_executable(test_factory_simple test_factory_simple.f90)
target_link_libraries(test_factory_simple
    PRIVATE 
        core
        infrastructure
        reconstructor
        flux
)

add_executable(test_domain_solution test_domain_solution.f90)
target_link_libraries(test_domain_solution
    PRIVATE 
        infrastructure
        core
)

add_executable(test_component_manager_physics test_component_manager_physics.f90)
target_link_libraries(test_component_manager_physics
    PRIVATE 
        manager
        infrastructure
        physics
        core
)
