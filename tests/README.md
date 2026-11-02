# OneFLOW GoogleTest tests

## 目录职责与来源

`tests/` 是 CMake/CTest 驱动的 GoogleTest 单元测试与 backend contract test 入口；它不替代既有的 `test/` 目录。`test/` 保留 OneFLOW 的 residual、算例和集成回归体系，`tests/` 面向小规模、确定性的 C++ 接口测试。

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
