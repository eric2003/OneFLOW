# src/core/CMakeLists.txt
# 创建核心库
add_library(core
    registry.f90
)

# 设置模块输出目录
target_include_directories(core PUBLIC
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# 安装规则（可选）
install(TARGETS core
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)