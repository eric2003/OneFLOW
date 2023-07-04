Eigen3
==================================

Eigen3
---------------------------------
#. `CMake+vcpkg编译简单Eigen3代码 <https://zhuanlan.zhihu.com/p/410353438/>`_
#. `Basic Operations in Eigen3 C++ <https://geophydog.cool/post/eigen3_operations/>`_



install
  
::

  .\vcpkg.exe install eigen3 
  
cmake  

::
  
  # this is heuristically generated, and may not be correct
  find_package(Eigen3 CONFIG REQUIRED)
  target_link_libraries(main PRIVATE Eigen3::Eigen)  
  
::
  
  cmake -DCMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" ../ 
  
::

  cmake -DCMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"  -DCMAKE_BUILD_TYPE=Release ../ 



