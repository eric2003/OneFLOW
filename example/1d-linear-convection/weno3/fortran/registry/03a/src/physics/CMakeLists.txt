# src/physics/CMakeLists.txt
message(STATUS "配置物理模块...")

# 创建物理模块库
add_library(physics STATIC
    physics_interface.f90
    equations/linear_convection.f90
    problems/linear_convection_problem.f90
)

# 链接依赖
target_link_libraries(physics PRIVATE base)

# 设置模块输出目录
set_target_properties(physics PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "物理模块配置完成")