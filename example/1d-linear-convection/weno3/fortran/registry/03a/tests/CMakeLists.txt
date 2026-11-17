# tests/CMakeLists.txt
message(STATUS "Configuring tests...")

#
#message(STATUS "CMAKE_Fortran_MODULE_DIRECTORY=${CMAKE_Fortran_MODULE_DIRECTORY}")
#
add_executable(test_minimal_simple test_minimal_simple.f90)
target_link_libraries(test_minimal_simple
    PRIVATE 
        core           # 必须链接core库
        infrastructure
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


add_executable(test_component_manager test_component_manager.f90)

target_link_libraries(test_component_manager
    PRIVATE 
        manager  # ← 链接到新的管理器库
        infrastructure
)
add_executable(test_basic_only test_basic_only.f90)
target_link_libraries(test_basic_only
    PRIVATE 
        infrastructure
        core
)

add_executable(test_physics_minimal test_physics_minimal.f90)
target_link_libraries(test_physics_minimal
    PRIVATE 
        physics
        base
)

add_executable(test_domain_solution test_domain_solution.f90)
target_link_libraries(test_domain_solution
    PRIVATE 
        infrastructure
        core
)

add_executable(test_config_physics test_config_physics.f90)
target_link_libraries(test_config_physics
    PRIVATE 
        infrastructure
        core
)

add_executable(test_solver_base test_solver_base.f90)
target_link_libraries(test_solver_base
    PRIVATE
        solver
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