# src/manager/CMakeLists.txt
message(STATUS "配置管理器模块...")

# 创建管理器库
add_library(manager STATIC
    component_manager.f90
    component_factory.f90
)

# 明确依赖关系：管理器依赖所有其他模块
target_link_libraries(manager
    PRIVATE 
        core
        infrastructure
        reconstructor
        flux
)

# 设置模块输出目录
set_target_properties(manager PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "管理器模块配置完成")