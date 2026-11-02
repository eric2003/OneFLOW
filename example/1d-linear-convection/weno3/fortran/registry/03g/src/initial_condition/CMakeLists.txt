# src/initial_condition/CMakeLists.txt
message(STATUS "配置初始条件模块...")

add_library(initial_condition STATIC
    ic_base.f90           # 基类
    step.f90              # 阶跃函数
    sine.f90              # 正弦波
    gaussian.f90          # 高斯脉冲
    factory.f90           # 工厂
)

target_link_libraries(initial_condition
    PRIVATE
        base
        infrastructure
        core                # 需要注册系统
)

set_target_properties(initial_condition PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "初始条件模块配置完成")