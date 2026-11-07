# OneFLOW GoogleTest tests

## 目录职责与来源

`tests/` 是 CMake/CTest 驱动的 GoogleTest 单元测试与 backend contract test 入口；它不替代既有的 `test/` 目录。`test/` 保留 OneFLOW 的 residual、算例和集成回归体系，`tests/` 面向小规模、确定性的 C++ 接口测试。JSONCpp 的解析、错误报告和序列化回环测试位于 `tests/io/`，使用内存流，不依赖当前工作目录中的临时 JSON 文件。

该目录随 GoogleTest 集成由提交 `8525ad1e`（2026-09-02）首次创建；提交 `60681aaa`（2026-09-02）加入 Euler backend contract test。当前包含原始 `hello_test` smoke test 和 CPU Euler backend contract test；Kunshan HIP contract test 由对应 port 的 CMake 选项按需构建。

默认行为是 CPU-only：根目录 CMake/CTest 不会申请或调用海光 DCU。Kunshan port 的 `ONEFLOW_1D_ENABLE_HIP` 和 `ONEFLOW_1D_ENABLE_GTEST` 默认均为 `OFF`，只有专用计算节点作业显式设置为 `ON` 时才构建并运行 HIP/DCU contract test。

mkdir build && cd build
cmake ..
cmake --build .
ctest --output-on-failure

# 查看是否有测试被注册
ctest -N -C Release

# 查看详细信息
ctest -N -C Release -V

PS D:\github\OneFLOW\build> ctest -N -C Release
Test project D:/github/OneFLOW/build
  Test #1: HelloTest.BasicAssertions
  Test #2: HelloTest.AnotherTest

Total Tests: 2

PS D:\github\OneFLOW\build> ctest --output-on-failure -C Release
Test project D:/github/OneFLOW/build
    Start 1: HelloTest.BasicAssertions
1/2 Test #1: HelloTest.BasicAssertions ........   Passed    0.10 sec
    Start 2: HelloTest.AnotherTest
2/2 Test #2: HelloTest.AnotherTest ............   Passed    0.03 sec

100% tests passed out of 2

Total Test time (real) =   0.28 sec

PS D:\github\OneFLOW\build> ctest -N -C Release -V
UpdateCTestConfiguration  from :D:/github/OneFLOW/build/DartConfiguration.tcl
Test project D:/github/OneFLOW/build
Constructing a list of tests
Done constructing a list of tests
Updating test list for fixtures
Added 0 tests to meet fixture requirements

1: Test command: "C:\Program Files\CMake\bin\cmake.exe" "-D" "TEST_EXECUTABLE=D:/github/OneFLOW/build/bin/Release/hello_test.exe" "-D" "TEST_EXECUTOR=" "-D" "TEST_FILTER=HelloTest.BasicAssertions" "-D" "TEST_XML_OUTPUT=" "-D" "TEST_EXTRA_ARGS=" "-P" "C:/Program Files/CMake/share/cmake-4.4/Modules/GoogleTest/LaunchTest.cmake"
1: Working Directory: D:/github/OneFLOW/build/tests
  Test #1: HelloTest.BasicAssertions

2: Test command: "C:\Program Files\CMake\bin\cmake.exe" "-D" "TEST_EXECUTABLE=D:/github/OneFLOW/build/bin/Release/hello_test.exe" "-D" "TEST_EXECUTOR=" "-D" "TEST_FILTER=HelloTest.AnotherTest" "-D" "TEST_XML_OUTPUT=" "-D" "TEST_EXTRA_ARGS=" "-P" "C:/Program Files/CMake/share/cmake-4.4/Modules/GoogleTest/LaunchTest.cmake"
2: Working Directory: D:/github/OneFLOW/build/tests
  Test #2: HelloTest.AnotherTest

Total Tests: 2

ctest -C Release -V --output-on-failure


## Production single-cell CFD smoke test

When CGNS support is enabled, CTest also registers
`GmshSingleCellRuntimeIntegration.ConvertsReadsComputesAndWritesFlowField`.
It creates an isolated temporary case, converts a one-cell 2D Gmsh QUAD_4 mesh
through a real CGNS interchange file into the native OneFLOW grid, then invokes
the production `OneFLOW r <case-directory>` command-line path for one solver
step. The test requires residual, restart-flow, and visualization outputs and
checks that the visualization contains exactly one cell.

This is a workflow smoke test, not a numerical verification case: it does not
compare the computed solution against an analytical solution or a CFD
benchmark. Its purpose is to expose mesh-read, single-cell/ghost-boundary,
solver-initialization, iteration, restart, and field-output integration defects.
All generated files are kept in the operating-system temporary directory.


## CMake infrastructure tests

The cmake/ subdirectory contains configure-level tests for the build system
itself. These are CTest tests, not GoogleTest cases and not production C++
tests.

CMakeSourceDiscovery creates a disposable fixture project below the build
tree and checks legacy and include/src layouts, nested source discovery,
header/include discovery, empty modules, preservation of sources under a
module directory named test, exclusion of generated build/CMakeFiles sources,
and source additions/removals across reconfiguration. It uses
project(... NONE), so this focused test does not require a C++ compiler or
third-party CFD libraries.

CMakeLegacyDiscoveryParity compares the new helper's source, header, and
include-directory sets with the existing ConstructSolutionDirTree function
against the real codes/ tree. This is the required compatibility gate before
migrating production callers.

Run only this test with:

.. code-block:: console

   ctest --test-dir build -R CMakeSourceDiscovery --output-on-failure

The fixture must remain isolated in the binary tree. It must not create or
modify files below the OneFLOW source tree.
