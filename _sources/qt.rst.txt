Qt
==================================

Qt is a cross-platform application development framework that allows developers to create high-performance, scalable, and portable applications. It provides a set of libraries and tools for various aspects such as GUI, networking, and databases. Qt also has an easy-to-use integrated development environment (IDE) called Qt Creator, which helps developers to develop and debug Qt applications more quickly.

#. `Qt's official website <https://www.qt.io/>`_
#. `Qt documentation <https://doc.qt.io/>`_
#. `QT Series Link Compilation <https://zhuanlan.zhihu.com/p/565066693/>`_
#. `Qt 6.3.1 C++ GUI Development Tutorial <https://zhuanlan.zhihu.com/p/565557087/>`_
#. `Qt6+windeployqt Series Link Compilation <https://zhuanlan.zhihu.com/p/566839520/>`_
#. `QMYSQL driver 6.6.1 <https://github.com/thecodemonkey86/qt_mysql_driver/>`_
#. `《Qt 5.9 C++开发指南》2021 完整版 <https://www.bilibili.com/video/BV1AX4y1w7Nt/>`_
#. `Qt6.3.1 C++ GUI开发教程（完整版） <https://www.bilibili.com/video/BV1G94y1Q7h6/>`_


Building from source (Qt6/CMake)
-----------------------------------
::

  set PATH=%PATH%;C:\Qt\Tools\CMake_64\bin;C:\Qt\Tools\Ninja
  C:
  cd C:\Qt\6.6.0\Src\qtbase\src\plugins\sqldrivers
  call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
  call C:\Qt\6.6.0\msvc2019_64\bin\qt-cmake.bat -G "Ninja Multi-Config" . -DMySQL_INCLUDE_DIR="E:\qt_creator\libs\libmysql\include" -DMySQL_LIBRARY="E:\qt_creator\libs\libmysql\lib\libmysql.lib" -DCMAKE_INSTALL_PREFIX="C:\Qt\6.6.0\msvc2019_64" -DCMAKE_CONFIGURATION_TYPES=Release;Debug
  ninja
  ninja install
  pause
  
Install Qt by MaintenanceTool.exe
--------------------------------------
::
  
  MaintenanceTool.exe


View QT version number
-----------------------
::

  ./qmake -v
  results:
  PS C:\local\Qt\6.6.1\msvc2019_64\bin> ./qmake -v
  QMake version 3.1
  Using Qt version 6.6.1 in C:/local/Qt/6.6.1/msvc2019_64/lib






