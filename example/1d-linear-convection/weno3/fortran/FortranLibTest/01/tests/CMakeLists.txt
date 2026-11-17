# 构建测试可执行文件（依赖mymath库）
add_executable(test_mymath test_mymath.f90)

# 关键：链接mymath库（自动继承mymath的公共包含目录，无需额外指定模块目录）
target_link_libraries(test_mymath PUBLIC mymath)


