Eigen3
==================================

Eigen3
---------------------------------
#. `CMake+vcpkg编译简单Eigen3代码 <https://zhuanlan.zhihu.com/p/410353438/>`_
#. `Basic Operations in Eigen3 C++ <https://geophydog.cool/post/eigen3_operations/>`_
#. `Introduction to Eigen C++ Matrix Library <https://aleksandarhaber.com/starting-with-eigen-c-matrix-library/>`_
#. `Eigen Tutorial 中文文档(c++版) <https://zhuanlan.zhihu.com/p/87613088/>`_
#. `C++ and Eigen Tutorials <https://www.youtube.com/playlist?list=PLO89phzZmnHjawqmeIbxXyIIZxhfgxut5/>`_
#. `Matrix Manipulations in C++ using Eigen Library <https://iamfaisalkhan.com/matrix-manipulations-using-eigen-cplusplus/>`_
#. `分布式并行计算笔记-MPI+openmp+Eigen <https://zhuanlan.zhihu.com/p/573503615/>`_
#. `上海交大超算平台用户手册 Eigen <https://docs.hpc.sjtu.edu.cn/app/compilers_and_languages/eigen.html>`_
#. `Kalman Filter: How to Implement in C++ with Eigen <https://codingcorner.org/kalman-filter-cpp-eigen-cmake/>`_
#. `Eigen 3.3.7 <https://apolo-docs.readthedocs.io/en/latest/software/scientific_libraries/eigen/eigen-3.3.7/index.html>`_
#. `Eigen 3.4.0 devdocs <https://devdocs.io/eigen3/>`_



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



