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

message(STATUS "源代码目录配置完成")