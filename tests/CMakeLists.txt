
# 包含 GoogleTest 模块 (只需在顶层包含一次，子目录即可直接使用)
include(GoogleTest)

# CMake infrastructure tests run in isolated temporary fixture projects.
add_subdirectory(cmake)

if ( NOT COMMAND gtest_discover_tests )
  message(FATAL_ERROR
    "GoogleTest.cmake did not define gtest_discover_tests. "
    "Check the CMake installation and ensure include(GoogleTest) "
    "is loaded before adding tests.")
endif()

# Optional external-tool integration tests.
# AUTO keeps the core test suite independent from optional executables.
# ON requires the currently supported external tool (Gmsh); OFF disables
# external-tool integration tests even when the tool is installed.
set(ONEFLOW_EXTERNAL_TOOL_TESTS "AUTO" CACHE STRING
    "External-tool tests: AUTO (detect), ON (require), or OFF (disable)")
set_property(CACHE ONEFLOW_EXTERNAL_TOOL_TESTS PROPERTY STRINGS AUTO ON OFF)
string(TOUPPER "${ONEFLOW_EXTERNAL_TOOL_TESTS}" _oneflow_external_tool_tests)

if (NOT _oneflow_external_tool_tests MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR
    "ONEFLOW_EXTERNAL_TOOL_TESTS must be AUTO, ON, or OFF; got "
    "'${ONEFLOW_EXTERNAL_TOOL_TESTS}'.")
endif()

find_program(ONEFLOW_GMSH_EXECUTABLE NAMES gmsh gmsh.exe)
set(ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS OFF)

if (_oneflow_external_tool_tests STREQUAL "ON")
  if (NOT ONEFLOW_GMSH_EXECUTABLE)
    message(FATAL_ERROR
      "ONEFLOW_EXTERNAL_TOOL_TESTS=ON requires the Gmsh executable on PATH.")
  endif()
  set(ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS ON)
elseif (_oneflow_external_tool_tests STREQUAL "AUTO")
  if (ONEFLOW_GMSH_EXECUTABLE)
    set(ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS ON)
  endif()
endif()

if (ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS)
  message(STATUS
    "External-tool tests: Gmsh integration enabled (${ONEFLOW_GMSH_EXECUTABLE})")
else()
  message(STATUS "External-tool tests: Gmsh integration disabled")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# 设置默认文件夹名称
set(CMAKE_FOLDER "Tests")

# L0 mesh metadata and validation tests with no solver or geometry-core dependency.
add_subdirectory(grid_summary)

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

# 4. database 模块测试 
add_subdirectory(database)

# 5. task 模块测试 
add_subdirectory(task)

# 6. register 模块测试 
add_subdirectory(register)

# 7. main / simulation entry & TaskRegistry (phase 0+1)
add_subdirectory(main)

# 8. zone
add_subdirectory(zone)

# 9. solver
add_subdirectory(solver)

# 10. io
add_subdirectory(io)

# 11. geometry
add_subdirectory(geometry)

# 12. basic utilities
add_subdirectory(basic)

# 13. accelerator runtime lifecycle
add_subdirectory(accel)

