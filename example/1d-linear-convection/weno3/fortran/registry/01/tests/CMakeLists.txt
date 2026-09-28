# tests/CMakeLists.txt
# 创建测试程序
add_executable(test_registry
    test_registry.f90
)

# 链接核心库
target_link_libraries(test_registry
    core
)

# 设置模块路径
target_include_directories(test_registry PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)