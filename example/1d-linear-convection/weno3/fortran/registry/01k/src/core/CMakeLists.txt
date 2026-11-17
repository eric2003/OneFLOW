# src/core/CMakeLists.txt
message(STATUS "配置核心模块...")

add_library(core STATIC
    registry.f90
    factory_interfaces.f90
)

target_include_directories(core PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# 安装规则
install(TARGETS core
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
install(FILES
    ${CMAKE_Fortran_MODULE_DIRECTORY}/registry_module.mod
    ${CMAKE_Fortran_MODULE_DIRECTORY}/factory_interfaces.mod
    DESTINATION include/fortran_cfd/core
)

message(STATUS "核心模块配置完成")