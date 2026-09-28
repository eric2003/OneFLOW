# 修改后的完整版本：
message(STATUS "配置源代码目录...")

# 设置包含目录
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# 按依赖顺序添加子目录
add_subdirectory(base)         # 1. 基础类型和精度
add_subdirectory(core)         # 2. 核心注册系统
add_subdirectory(infrastructure) # 3. 基础设施（配置、网格、域、解）
add_subdirectory(physics)      # 4. 物理模块（方程和问题）
add_subdirectory(numerics)     # 5. 数值方法（重构器、通量、时间积分）
add_subdirectory(manager)      # 6. 组件管理器和工厂
add_subdirectory(solver)       # 7. 求解器

# ==================== 新增：边界条件和初始条件模块 ====================
message(STATUS "配置边界条件模块...")
add_subdirectory(boundary)

message(STATUS "配置初始条件模块...")
add_subdirectory(initial_condition)

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

message(STATUS "集成求解器库配置完成")