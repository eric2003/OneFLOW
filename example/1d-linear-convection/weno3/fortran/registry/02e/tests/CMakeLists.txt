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

# 更新测试链接
#add_executable(test_component_manager test_component_manager.f90)
#
#target_link_libraries(test_component_manager
#    PRIVATE 
#        manager  # ← 链接到新的管理器库
#        infrastructure
#)
#
#add_executable(test_architecture test_architecture.f90)
#target_link_libraries(test_architecture
#    PRIVATE 
#        core
#        infrastructure
#        manager
#)

#add_executable(test_solver_framework test_solver_framework.f90)
#target_link_libraries(test_solver_framework
#    PRIVATE 
#        core
#        infrastructure
#        manager
#)
#
#add_executable(test_cfd_architecture test_cfd_architecture.f90)
#target_link_libraries(test_cfd_architecture
#    PRIVATE
#        base
#        core
#        infrastructure
#)


add_executable(test_basic_only test_basic_only.f90)
target_link_libraries(test_basic_only
    PRIVATE 
        infrastructure
        core
)