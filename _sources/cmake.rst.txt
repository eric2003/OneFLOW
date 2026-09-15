CMake
==================================

CMake is a cross-platform build system generator. It allows you to define a build process in a high-level language and generates native build files such as Makefiles or Visual Studio projects for different platforms and compilers. This allows you to build your software in a consistent and reproducible way across multiple platforms and development environments.

CMake uses a simple, declarative syntax that is easy to learn and understand. You define your project's source files, libraries, and executable targets in a CMakeLists.txt file, and CMake takes care of the rest. You can also use CMake to configure your project with different build options, such as enabling or disabling certain features or selecting different build configurations.

CMake supports a wide variety of platforms and compilers, including Windows, Linux, macOS, and many embedded platforms. It can generate build files for a variety of build systems, including Make, Ninja, and Visual Studio. CMake can also integrate with popular IDEs such as Visual Studio, Eclipse, and Xcode, making it easy to use with your favorite development tools.

Overall, CMake provides a powerful and flexible build system that can simplify your build process and make your software more portable and maintainable.

Generate a Project Buildsystem
::

 cmake [<options>] <path-to-source | path-to-existing-build>
 cmake [<options>] -S <path-to-source> -B <path-to-build>

Build a Project
::

 cmake --build <dir> [<options>] [-- <build-tool-options>]

Install a Project
::

 cmake --install <dir> [<options>]

Open a Project
::

 cmake --open <dir>

Run a Script
::

 cmake [-D <var>=<value>]... -P <cmake-script-file>

Run a Command-Line Tool
::

 cmake -E <command> [<options>]

Run the Find-Package Tool
::

 cmake --find-package [<options>]

Run a Workflow Preset
::

 cmake --workflow [<options>]

View Help
::

 cmake --help[-<topic>]
 
CMAKE_PREFIX_PATH
::

  cmake .. -D CMAKE_PREFIX_PATH:PATH="c:/dev/googletest;c:/dev/cgns/"

#. `CMake official documentation <https://cmake.org/documentation/>`_
#. `CMake official tutorial <https://cmake.org/cmake/help/latest/guide/tutorial/>`_
#. `CMake official GitHub repository <https://github.com/Kitware/CMake/>`_
#. `Compilation of links related to CMake from beginner to expert series <https://zhuanlan.zhihu.com/p/393316878/>`_
#. `CMake hands-on workshop <https://enccs.github.io/cmake-workshop/>`_
#. `import CMake, CMake and C++20 Modules - Bill Hoffman - CppCon 2022 <https://www.youtube.com/watch?v=5X803cXe02Y/>`_
#. `import CMake: // 2023 State of C++20 modules in CMake - Bill Hoffman - CppNow 2023 <https://www.youtube.com/watch?v=c563KgO-uf4/>`_
#. `So, You Want to Use C++ Modules … Cross-Platform? - Daniela Engert - C++ on Sea 2023 <https://www.youtube.com/watch?v=DJTEUFRslbI/>`_
#. `C++ Modules in 2023 <https://www.youtube.com/watch?v=vAjEkIy43-c/>`_
#. `C++ Features You Might Not Know - Jonathan Müller - C++ on Sea 2023 <https://www.youtube.com/watch?v=zGWj7Qo_POY/>`_

Set Start project(Visual Studio)
::

  if (MSVC)
    set_property ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      PROPERTY 
         VS_STARTUP_PROJECT ${PROJECT_NAME}
    )
  endif () 

CMake+Qt Example
::

 cmake .. -D CMAKE_PREFIX_PATH:PATH=C:/local/Qt/Qt6.4.0/6.4.0/msvc2019_64/
 cmake --build .
 C:/local/Qt/Qt6.4.0/6.4.0/msvc2019_64/bin/windeployqt.exe .\Debug\testprj.exe
 
CMake+HDF5 Example(powershell)
::

 $env:HDF5_DIR = "c:/dev/HDF_Group/HDF5/1.14.0/cmake/"
 $env:path += ";C:/dev/HDF_Group/HDF5/1.14.0/bin/"
 cmake ..
 cmake --build .


Debug
::

  cmake -DCMAKE_BUILD_TYPE=Debug ..
  
Debug
::

  cmake --build . --parallel 4 --config Debug  
  
Release
::

  cmake -DCMAKE_BUILD_TYPE=Release ..  
  
Release
::

  cmake --build . --parallel 4 --config Release
  
Install
::

  cmake --install . --prefix "c:/dev/myprj/"  
  or
  cmake --install . --prefix "c:/dev/myprj/"  --config Debug
  
Display all current environment variables
::
 
  cmake -E environment  

$ENV
::
  
  message ( STATUS "ENV{OS} = $ENV{OS}" )  
  
Some system variables
::

  ProgramData                    C:\ProgramData
  ProgramFiles                   C:\Program Files
  ProgramFiles(x86)              C:\Program Files (x86)
  ProgramW6432                   C:\Program Files  
  
Display ProgramFiles(x86) variable
::
  
  message ( STATUS "ENV{ProgramFiles\(x86\)} = $ENV{ProgramFiles\(x86\)}" )
  
Powershell Remove
::

   Remove-Item * -Recurse -Force 
   
Save CMake output to file
::

  cmake .. >> output_file.txt 2>&1
  or 
  cmake --version | Set-Content -Path AAA.txt
  or 
  cmake --version | Out-File -FilePath C:\TestDir\AliasNames.txt

trace
::

  cmake --trace .. >> output.txt 2>&1
  
default trace
::

  CMakeDetermineSystem
  configure_file(${CMAKE_ROOT}/Modules/
  CMakeSystem.cmake.in ${CMAKE_PLATFORM_INFO_DIR}/CMakeSystem.cmake @ONLY )
  CMakeSystemSpecificInitialize
  CMakeDetermineCCompiler
  CMakeDetermineCompiler
  CMakeDetermineCompilerId
  CMakeCompilerIdDetection
  ARMCC-DetermineCompiler
  ARMClang-DetermineCompiler
  AppleClang-DetermineCompiler
  Clang-DetermineCompilerInternal
  Borland-DetermineCompiler
  Bruce-C-DetermineCompiler
  Clang-DetermineCompilerInternal
  Compaq-C-DetermineCompiler
  Cray-DetermineCompiler
  CrayClang-DetermineCompiler
  Embarcadero-DetermineCompiler
  Fujitsu-DetermineCompiler
  FujitsuClang-DetermineCompiler
  GHS-DetermineCompiler
  GNU-C-DetermineCompiler
  HP-C-DetermineCompiler
  IAR-DetermineCompiler
  IBMClang-C-DetermineCompiler
  Intel-DetermineCompiler
  IntelLLVM-DetermineCompiler
  NVIDIA-DetermineCompiler
  OrangeC-DetermineCompiler
  PathScale-DetermineCompiler
  SDCC-C-DetermineCompiler
  SunPro-C-DetermineCompiler
  TI-DetermineCompiler
  Tasking-DetermineCompiler
  VisualAge-C-DetermineCompiler
  IBMCPP-C-DetermineVersionInternal
  Watcom-DetermineCompiler
  XL-C-DetermineCompiler
  XLClang-C-DetermineCompiler
  zOS-C-DetermineCompiler
  CMakeFindBinUtils
  D:/work/modern_cmake_work/ModernCMake/codes/cmake/trace/01/build/CMakeFiles/3.28.0-rc4/CMakeCCompiler.cmake
  CMakeDetermineCXXCompiler
  D:/work/modern_cmake_work/ModernCMake/codes/cmake/trace/01/build/CMakeFiles/3.28.0-rc4/CMakeCXXCompiler.cmake
  CMakeSystemSpecificInformation
  CMakeGenericSystem
  CMakeInitializeConfigs
  Windows
  WindowsPaths
  CMakeCInformation
  CMakeLanguageInformation
  MSVC-C
  MSVC
  CMakeCommonCompilerMacros
  CMakeCInformation
  Windows-MSVC-C
  Windows-MSVC
  CMakeDetermineRCCompiler
  CMakeRCCompiler
  D:/work/modern_cmake_work/ModernCMake/codes/cmake/trace/01/build/CMakeFiles/3.28.0-rc4/CMakeRCCompiler.cmake
  CMakeRCInformation
  CMakeTestRCCompiler
  CMakeCommonLanguageInclude
  CMakeTestCCompiler
  CMakeTestCompilerCommon
  CMakeDetermineCompilerABI
  CMakeParseImplicitIncludeInfo
  CMakeParseLibraryArchitecture
  D:/work/modern_cmake_work/ModernCMake/codes/cmake/trace/01/build/CMakeFiles/CMakeScratch/TryCompile-6now2o/CMakeLists.txt
  Windows-Initialize
  CMakeDetermineCompileFeatures
  FeatureTesting
  
Find CUDA IMPORTED_TARGETS
::

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
  
  find_package ( CUDAToolkit )
  
  get_directory_property( my_import_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} IMPORTED_TARGETS )
  
  message( STATUS "my_import_targets=${my_import_targets}" ) 

results
::

  my_import_targets=CUDA::toolkit;CUDA::cuda_driver;CUDA::cudart;CUDA::cudart_static;CUDA::cudart_static_deps;
  CUDA::cublasLt;CUDA::cublas;CUDA::cufft;CUDA::curand;CUDA::cusparse;CUDA::nppc;CUDA::nvjpeg;CUDA::cufftw;
  CUDA::cusolver;CUDA::nppial;CUDA::nppicc;CUDA::nppidei;CUDA::nppif;CUDA::nppig;CUDA::nppim;CUDA::nppist;
  CUDA::nppitc;CUDA::npps;CUDA::nppisu;CUDA::cupti;CUDA::nvrtc;CUDA::nvml;CUDA::nvToolsExt;CUDA::OpenCL

vcpkg
::

  cmake -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" ..
  
  
CMake --help-property-list
-------------------------------
#. `CMake --help-property-list <https://zhuanlan.zhihu.com/p/536789898/>`_
  
generator_expressions
::

  add_custom_target ( print 
      ${CMAKE_COMMAND} -E 
      echo 
      BUILD_INTERFACE:a = $<BUILD_INTERFACE:a> &
      echo 
      INSTALL_INTERFACE:a = $<INSTALL_INTERFACE:a>
  )
  
target print
::

  cmake --build . --target print
  cmake --build . --config Debug --target print
  cmake --build . --config Release --target print
  
print_variables
::

  include(CMakePrintHelpers)
  cmake_print_variables(CMAKE_C_COMPILER CMAKE_MAJOR_VERSION DOES_NOT_EXIST)  
  
print_properties
::  
  
  include(CMakePrintHelpers)
  cmake_print_properties(
    TARGETS ${PROJECT_NAME}
    PROPERTIES LINK_LIBRARIES
  )
  
  cmake_print_properties(
    TARGETS ${PROJECT_NAME}
    PROPERTIES COMPILE_FEATURES
  )  
  
print_properties function
::

  function(print_varname inputVarName )
      message ( STATUS "inputVarName    = ${inputVarName}" )
      get_property( variableNames DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VARIABLES)
      list (SORT variableNames )
      foreach ( varName ${variableNames})
          string( TOUPPER ${varName} upperVarName )
          string ( FIND ${upperVarName} ${inputVarName} myloc )
          if ( myloc GREATER_EQUAL 0 )
              message( STATUS "${varName}" )
          endif () 
      endforeach()
  endfunction()
  
  function(print_varname_value inputVarName )
      message ( STATUS "inputVarName    = ${inputVarName}" )
      get_property( variableNames DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VARIABLES)
      list (SORT variableNames )
      foreach ( varName ${variableNames})
          string( TOUPPER ${varName} upperVarName )
          string ( FIND ${upperVarName} ${inputVarName} myloc )
          if ( myloc GREATER_EQUAL 0 )
              message( STATUS "${varName}=${${varName}}" )
          endif () 
      endforeach()
  endfunction()  
  
  print_varname( "INCLUDE" )
  print_varname( "LIB" )
  print_varname( "EIGEN" )  
  
print_target_properties:
::

  function(my_print_target_properties targetname )
      include(CMakePrintHelpers)
      message ( STATUS "targetname = ${targetname}" )
      set( props ${ARGN} )
      foreach ( prop IN LISTS props )
          cmake_print_properties(
            TARGETS ${targetname}
            PROPERTIES ${prop}
          )
      endforeach()
  endfunction()  

some property
::
  
  COMPILE_DEFINITIONS
  COMPILE_DEFINITIONS_<CONFIG>
  COMPILE_FEATURES
  COMPILE_FLAGS
  COMPILE_OPTIONS
  DEFINITIONS
  ENABLED_FEATURES
  ENABLED_LANGUAGES
  ENABLE_EXPORTS
  HEADER_DIRS
  IMPORTED
  
::

  function(print_arguments)
    foreach(arg ${ARGN})
      message("Argument: ${arg}")
    endforeach()
  endfunction()
  
  print_arguments("Hello" "World" "CMake")  
  
::

  function(print_var var)
    message( STATUS "ARGC: ${ARGC}" )
    message( STATUS "ARGV: ${ARGV}" )
    message( STATUS "ARGN: ${ARGN}" )
    foreach(arg IN LISTS ${ARGN})
      message("Argument: ${arg}")
    endforeach()
  endfunction()
  
  #print_var() #error
  print_var(a)
  print_var(a b)
  print_var(a b c)
  
results:  
::

  -- ARGC: 1
  -- ARGV: a
  -- ARGN:
  -- ARGC: 2
  -- ARGV: a;b
  -- ARGN: b
  -- ARGC: 3
  -- ARGV: a;b;c
  -- ARGN: b;c  
  
  
::  

  function(echo_target tgt)
      if(NOT TARGET ${tgt})
          message("There is no target named '${tgt}'")
          return()
      endif()
  
      set(props
          DEBUG_OUTPUT_NAME
          DEBUG_POSTFIX
          RELEASE_OUTPUT_NAME
          ...
          LINK_SEARCH_START_STATIC
          LOCATION
          LOCATION_DEBUG
          ...
          WIN32_EXECUTABLE
          XCODE_ATTRIBUTE_WHATEVER
      )
      message(STATUS "======================== ${tgt} ========================")
  
      # Push the current (NEW) CMake policy onto the stack, and apply the OLD policy.
      cmake_policy(PUSH)
      cmake_policy(SET CMP0026 OLD)
  
      foreach(p ${props})
          # v for value, d for defined, s for set
          get_property(v TARGET ${tgt} PROPERTY ${p})
          get_property(d TARGET ${tgt} PROPERTY ${p} DEFINED)
          get_property(s TARGET ${tgt} PROPERTY ${p} SET)
          # only produce output for values that are set
          if(s)
              message(STATUS "tgt='${tgt}' p='${p}'")
              message(STATUS "  value='${v}'")
              message(STATUS "  defined='${d}'")
              message(STATUS "  set='${s}'")
              message(STATUS "")
          endif()
      endforeach()
  
      # Pop the previous policy from the stack to re-apply the NEW behavior.
      cmake_policy(POP)
  
      message(STATUS "")
  endfunction()
  
configure_file
::

  option(FOO_ENABLE "Enable Foo" ON)
  if(FOO_ENABLE)
    set(FOO_STRING "foo")
  endif()
  configure_file(foo.h.in foo.h @ONLY)

foo.h.in
::

  #cmakedefine FOO_ENABLE
  #cmakedefine FOO_STRING "@FOO_STRING@"
  
foo.h
::

  #define FOO_ENABLE
  #define FOO_STRING "foo"
  
configure_file example 2
::

  set(A_ENABLE ON)
  set(B_ENABLE ON)
  set(C_ENABLE OFF)
  set(FOO_STRING "foo")
  configure_file(foo.h.in foo.h @ONLY)

foo.h.in
::

  #cmakedefine A_ENABLE
  #cmakedefine01 B_ENABLE
  #cmakedefine01 C_ENABLE
  #cmakedefine FOO_STRING "@FOO_STRING@"
  
foo.h
::

  #define A_ENABLE
  #define B_ENABLE 1
  #define C_ENABLE 0
  #define FOO_STRING "foo"

Foo lib
::

  d:\work\github_work\Foo\
  
CMake+importing-exporting
---------------------------------------
#. `CMake+importing-exporting系列链接整理 <https://zhuanlan.zhihu.com/p/533480192/>`_

CMAKE_PREFIX_PATH
::

  cmake -D CMAKE_PREFIX_PATH="C:/local/MyInstall/MathFunctions" ..
  
Import target "Foo::Foo" for configuration "Debug"
::

  set_property(TARGET Foo::Foo APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)  
  
DCMAKE_VERBOSE_MAKEFILE
::

  cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ..
  
CMake FetchContent vs. ExternalProject
---------------------------------------
#. `CMake FetchContent vs. ExternalProject <https://www.scivision.dev/cmake-fetchcontent-vs-external-project/>`_

add_test
::

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
   
  add_executable( testprj )
  
  target_sources( testprj 
      PRIVATE 
        main.cpp
  )
  
  # Enable testing
  enable_testing()
  
  # Add a test
  add_test(
      NAME MyTest
      COMMAND testprj
  )

ctest
::

  cmake ..
  cmake --build .
  ctest -C Debug
  
add_test example 1
main.cpp
::

  #include <iostream>
  #include <cmath>
  #include <string>
  
  int main( int argc, char ** argv )
  {
      if ( argc < 2 ) {
        // report version
        std::cout << argv[0] << std::endl;
        std::cout << "Usage: " << argv[0] << " number" << std::endl;
        return 1;
      }
      
      // convert input to double
      const double inputValue = std::stod(argv[1]);
      
      const double outputValue = std::sqrt(inputValue);
      
      std::cout << "The square root of " << inputValue << " is " << outputValue
                << std::endl;
      return 0;
  }
  
CMakeList.txt
::

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
   
  add_executable( testprj )
  
  target_sources( testprj 
      PRIVATE 
        main.cpp
  )
  
  # Enable testing
  enable_testing()
  
  # does the usage message work?
  add_test(NAME Usage COMMAND testprj)
  # Add a test
  set_tests_properties(Usage
    PROPERTIES PASS_REGULAR_EXPRESSION "Usage:.*number"
  )
  
  
ctest
::

  cmake ..
  cmake --build .
  ctest -C Debug
  results:
  PS D:\work\modern_cmake_work\ModernCMake\codes\cmake\add_test\04\build> ctest -C Debug
  Test project D:/work/modern_cmake_work/ModernCMake/codes/cmake/add_test/04/build
      Start 1: Usage
  1/1 Test #1: Usage ............................   Passed    0.04 sec
  
  100% tests passed, 0 tests failed out of 1
  
  Total Test time (real) =   0.04 sec
  
::

  当使用ctest命令时，您可以使用这些参数来控制测试的行为。以下是一些具体的例子：
  
  使用EXCLUDE_REGEX参数排除特定的测试用例：
  
  ctest -E "slow_tests"  # 排除名称中包含"slow_tests"的测试用例
  使用INCLUDE_REGEX参数只运行特定的测试用例：
  
  ctest -R "test_[0-9]"  # 只运行名称以"test_"开头并且后面跟着一个数字的测试用例
  使用TIMEOUT参数设置测试用例的超时时间：
  
  ctest --timeout 300  # 设置测试用例的超时时间为300秒
  使用LABELS参数标记测试用例：
  
  ctest -L "unit_tests"  # 只运行标记为"unit_tests"的测试用例
  这些参数可以单独使用，也可以组合使用来满足具体的测试需求。您可以根据项目的实际情况选择适合的参数来控制测试的行为。  
  
Unit testing with CTest
---------------------------------------
#. `Unit testing with CTest <https://bertvandenbroucke.netlify.app/2019/12/12/unit-testing-with-ctest/>`_
#. `CMake的PASS_REGULAR_EXPRESSION如何匹配多行输出？ <https://cloud.tencent.com/developer/ask/sof/113372642/>`_
#. `Automated Testing with CMake, CTest and CDash <https://www.youtube.com/watch?v=YlIqlVVJWuo/>`_
#. `Do you even test? (your code with CMake) <https://www.youtube.com/watch?v=pxJoVRfpRPE/>`_
#. `Google Test and Mock Platform, Complete Tutorial. Part 1: Google Test <https://www.youtube.com/watch?v=JJqRlSTQlh4/>`_


googletest imported targets
::

  find_package(GTest REQUIRED)
  
  get_directory_property( my_import_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} IMPORTED_TARGETS )
  
  message( STATUS "my_import_targets=${my_import_targets}" )
  
  results:
  -- Found GTest: C:/dev/googletest/lib/cmake/GTest/GTestConfig.cmake (found version "1.14.0")
  -- my_import_targets=GTest::gtest;GTest::gtest_main;GTest::gmock;GTest::gmock_main;GTest::GTest;GTest::Main
  
  
declare our test, by specifying which command to run: 
::

  add_test(
    NAME cpp_test
    COMMAND $<TARGET_FILE:cpp_test>
  )
  
  
.. code-block:: cmake

   add_test(
     NAME testprj
     COMMAND $<TARGET_FILE:testprj>
   )  
  
cmake --build . --target print results:
::

  PS D:\work\modern_cmake_work\ModernCMake\codes\cmake\add_test\06\build> cmake --build . --target print
  适用于 .NET Framework MSBuild 版本 17.8.3+195e7f5a3
  
    testprj.vcxproj -> D:\work\modern_cmake_work\ModernCMake\codes\cmake\add_test\06\build\Debug\testprj.exe
    1>
    TARGET_FILE:testprj = D:/work/modern_cmake_work/ModernCMake/codes/cmake/add_test/06/build/Debug/testprj.exe  

导入库文件
使用add_library命令，通过指定IMPORTED选项表明这是一个导入的库文件，通过设置其属性指明其路径：

.. code-block:: cmake

  add_library(math STATIC IMPORTED)
  set_property(TARGET math PROPERTY
               IMPORTED_LOCATION "./lib/libmath.a")

对于库文件的路径，也可以使用find_library命令来查找，比如在lib目录下查找math的Realse和Debug版本：

.. code-block:: cmake

  find_library(LIB_MATH_DEBUG mathd HINTS "./lib")
  find_library(LIB_MATH_RELEASE math HINTS "./lib")
  
对于不同的编译类型，可以通过IMPORTED_LOCATION_<CONFIG>来指明不同编译类型对应的库文件路径：

.. code-block:: cmake

  add_library(math STATIC IMPORTED GLOBAL)
  set_target_properties(math PROPERTIES
    IMPORTED_LOCATION "${LIB_MATH_RELEASE}"
    IMPORTED_LOCATION_DEBUG "${LIB_MATH_DEBUG}"
    IMPORTED_CONFIGURATIONS "RELEASE;DEBUG"
  )
  
  
导入成功以后，就可以将该库链接到其他目标上，但是导入的目标不可以被install。
这里以导入静态库为例，导入动态库或其他类型也是类似的操作，只需要将文件类型STATIC修改成对应的文件类型即可。  

导入可执行文件
这个不是那么常用，为了文章完整性，顺便提一下。是和导入库文件类似的：

.. code-block:: cmake

  add_executable(demo IMPORTED)
  set_property(TARGET demo PROPERTY
               IMPORTED_LOCATION "./bin/demo")

bigobj problem

.. code-block:: cmake

  if ( MSVC )
    target_compile_options( testprj 
      PRIVATE
        /bigobj
    )
  endif()
  
string(REGEX REPLACE 

.. code-block:: cmake

  set(myString "Hello, World!")
  string(REGEX REPLACE "Hello" "Hi" myOutString ${myString})
  message ( STATUS "myString = ${myString}" )
  message ( STATUS "myOutString = ${myOutString}" )
  
results:
::

  -- myString = Hello, World!
  -- myOutString = Hi, World!  
  
  
在CMake中，使用message(TRACE ...)命令可以输出跟踪级别的消息。这些消息通常用于调试目的，但默认情况下，CMake可能不会显示跟踪级别的消息。

要查看跟踪级别的消息，可以在运行CMake时设置CMAKE_MESSAGE_LOG_LEVEL变量为TRACE。例如，在命令行中运行以下命令：  

::

  cmake -DCMAKE_MESSAGE_LOG_LEVEL=TRACE /path/to/your/source
  
message的高级使用-指定日志级别 message([<mode>] "message")
它的级别有--log-level = <ERROR|WARNING|NOTICE|STATUS|VERBOSE|DEBUG|TRACE>  

cmake -E
--------------
::

  Usage: C:\Program Files\CMake\bin\cmake.exe -E <command> [arguments...]
  Available commands:
    capabilities              - Report capabilities built into cmake in JSON format
    cat [--] <files>...       - concat the files and print them to the standard output
    chdir dir cmd [args...]   - run command in a given directory
    compare_files [--ignore-eol] file1 file2
                                - check if file1 is same as file2
    copy <file>... destination  - copy files to destination (either file or directory)
    copy_directory <dir>... destination   - copy content of <dir>... directories to 'destination' directory
    copy_directory_if_different <dir>... destination   - copy changed content of <dir>... directories to 'destination' directory
    copy_if_different <file>... destination  - copy files if it has changed
    echo [<string>...]        - displays arguments as text
    echo_append [<string>...] - displays arguments as text but no new line
    env [--unset=NAME ...] [NAME=VALUE ...] [--] <command> [<arg>...]
                              - run command in a modified environment
    environment               - display the current environment
    make_directory <dir>...   - create parent and <dir> directories
    md5sum <file>...          - create MD5 checksum of files
    sha1sum <file>...         - create SHA1 checksum of files
    sha224sum <file>...       - create SHA224 checksum of files
    sha256sum <file>...       - create SHA256 checksum of files
    sha384sum <file>...       - create SHA384 checksum of files
    sha512sum <file>...       - create SHA512 checksum of files
    remove [-f] <file>...     - remove the file(s), use -f to force it (deprecated: use rm instead)
    remove_directory <dir>... - remove directories and their contents (deprecated: use rm instead)
    rename oldname newname    - rename a file or directory (on one volume)
    rm [-rRf] [--] <file/dir>... - remove files or directories, use -f to force it, r or R to remove directories and their contents recursively
    sleep <number>...         - sleep for given number of seconds
    tar [cxt][vf][zjJ] file.tar [file/dir1 file/dir2 ...]
                              - create or extract a tar or zip archive
    time command [args...]    - run command and display elapsed time
    touch <file>...           - touch a <file>.
    touch_nocreate <file>...  - touch a <file> but do not create it.
    create_symlink old new    - create a symbolic link new -> old
    create_hardlink old new   - create a hard link new -> old
    true                      - do nothing with an exit code of 0
    false                     - do nothing with an exit code of 1
  Available on Windows only:
    delete_regv key           - delete registry value
    env_vs8_wince sdkname     - displays a batch file which sets the environment for the provided Windows CE SDK installed in VS2005
    env_vs9_wince sdkname     - displays a batch file which sets the environment for the provided Windows CE SDK installed in VS2008
    write_regv key value      - write registry value

Fortran CMake
::

  cmake_minimum_required(VERSION 3.28)
  
  project ( testprj )
  
  set ( PRJ_COMPILE_FEATURES )
  set ( PRJ_COMPILE_DEFINITIONS )
  
  enable_language(Fortran)
  
  if (MSVC)
    set_property ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      PROPERTY 
         VS_STARTUP_PROJECT ${PROJECT_NAME}
    )
  endif ()  
  
  add_executable( ${PROJECT_NAME}
      main.f90
  )
  
  target_compile_features ( ${PROJECT_NAME} 
  	PRIVATE 
  		${PRJ_COMPILE_FEATURES}
  )
  
  target_compile_definitions ( ${PROJECT_NAME}
  	PRIVATE
  	   ${PRJ_COMPILE_DEFINITIONS} 
  )
