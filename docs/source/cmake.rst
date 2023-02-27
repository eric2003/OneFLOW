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

#. `CMake official documentation <https://cmake.org/documentation/>`_
#. `CMake official tutorial <https://cmake.org/cmake/help/latest/guide/tutorial/>`_
#. `CMake official GitHub repository <https://github.com/Kitware/CMake/>`_
#. `Compilation of links related to CMake from beginner to expert series <https://zhuanlan.zhihu.com/p/393316878/>`_



