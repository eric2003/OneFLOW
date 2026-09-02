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
