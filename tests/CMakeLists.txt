
# 包含 GoogleTest 模块 (只需在顶层包含一次，子目录即可直接使用)
include(GoogleTest)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# 设置默认文件夹名称
set(CMAKE_FOLDER "Tests")

# ============================================================================
# 模块级测试路由 (使用 add_subdirectory 引入各个子目录的 CMakeLists.txt)
# 这样顶层文件永远保持简洁，新增模块只需加一行 add_subdirectory
# ============================================================================

# 1. Hello 模块测试
add_subdirectory(hello)

# 2. ADT 模块测试
add_subdirectory(adt)

# 3. Euler 模块测试 
add_subdirectory(euler)


