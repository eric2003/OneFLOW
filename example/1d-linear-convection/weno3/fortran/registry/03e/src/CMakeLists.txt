# ==================== src\CMakeLists.txt ====================

# src/CMakeLists.txt
message(STATUS "配置源代码目录...")

# 设置包含目录
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# 添加子目录
add_subdirectory(base)
add_subdirectory(core)
add_subdirectory(infrastructure)
add_subdirectory(numerics)
add_subdirectory(physics)  # ← 新增物理模块目录
add_subdirectory(manager)
add_subdirectory(solver)

# ==================== 新增：结果模块 ====================
message(STATUS "配置结果模块...")

add_library(results STATIC
    results.f90  # 新增文件
)

target_link_libraries(results
    PRIVATE
        base
        infrastructure
        solver
)

set_target_properties(results PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(STATUS "结果模块配置完成")
# ==================== 结束新增 ====================

message(STATUS "源代码目录配置完成")