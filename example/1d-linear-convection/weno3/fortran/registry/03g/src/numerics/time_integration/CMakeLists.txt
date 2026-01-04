# src/numerics/time_integration/CMakeLists.txt
message(STATUS "配置时间积分模块...")

add_library(time_integration STATIC
    base.f90
    rk1.f90
    rk2.f90
    rk3.f90
    factory.f90
)

target_link_libraries(time_integration
    PRIVATE
        infrastructure
        base
        manager
        solver
)

set_target_properties(time_integration PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "时间积分模块配置完成")