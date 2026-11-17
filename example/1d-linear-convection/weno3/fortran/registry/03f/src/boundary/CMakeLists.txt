# src/boundary/CMakeLists.txt
message(STATUS "配置边界条件模块...")

add_library(boundary STATIC
    boundary_base.f90        # 基类
    periodic.f90             # 周期性边界
    dirichlet.f90           # Dirichlet边界
    neumann.f90             # Neumann边界
    factory.f90             # 工厂
)

target_link_libraries(boundary
    PRIVATE
        base
        infrastructure
        core                 # 需要注册系统
)

set_target_properties(boundary PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "边界条件模块配置完成")