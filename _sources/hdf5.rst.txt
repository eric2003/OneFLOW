HDF5
==================================

#. `CGNS and HDF5 compiling and linking related article links compilation <https://zhuanlan.zhihu.com/p/452874893/>`_

powershell add environmental path
::

  $env:path += ";C:/dev/HDF_Group/HDF5/1.14.2/bin/"
  
::
  
  https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/bin/windows/hdf5-1.14.3-Std-win10_64-vs17.zip
  
::

  $hdf5_major = 1
  $hdf5_minor = 14
  $hdf5_patch = 3
  $hdf5_dir1 = "hdf5-$hdf5_major.$hdf5_minor"
  $hdf5_dir2 = "$hdf5_dir1.$hdf5_patch"
  $hdf5_url = "https://support.hdfgroup.org/ftp/HDF5/releases/$hdf5_dir1/$hdf5_dir2/bin/windows/"  
  
HDF5 Source Code
::

  https://www.hdfgroup.org/downloads/hdf5/source-code/  
  
  
CMake message

.. code-block:: cmake

  message ( STATUS "HDF_CONFIG_DIR = ${HDF_CONFIG_DIR}" )
  message ( STATUS "HDF_RESOURCES_DIR = ${HDF_RESOURCES_DIR}" )
  message ( STATUS "HDF5_SOURCE_DIR = ${HDF5_SOURCE_DIR}" )
  message ( STATUS "HDF5_SRC_DIR = ${HDF5_SRC_DIR}" )
  message ( STATUS "HDF5_TEST_SRC_DIR = ${HDF5_TEST_SRC_DIR}" )
  message ( STATUS "HDF5_TEST_PAR_DIR = ${HDF5_TEST_PAR_DIR}" )
  message ( STATUS "HDF5_TEST_API_SRC_DIR = ${HDF5_TEST_API_SRC_DIR}" )
  message ( STATUS "HDF5_TEST_API_PAR_SRC_DIR = ${HDF5_TEST_API_PAR_SRC_DIR}" )  
  
HDF5 CMake Info
::

  -- CONFIG_DATE = 2023-12-18
  -- HDF5_USE_FOLDERS = ON
  -- HDF_CONFIG_DIR = D:/work/hdf5_work/hdf5-1.14.3/config
  -- HDF_RESOURCES_DIR = D:/work/hdf5_work/hdf5-1.14.3/config/cmake
  -- HDF5_SOURCE_DIR = D:/work/hdf5_work/hdf5-1.14.3
  -- HDF5_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/src
  -- HDF5_TEST_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/test
  -- HDF5_TEST_PAR_DIR = D:/work/hdf5_work/hdf5-1.14.3/testpar
  -- HDF5_TEST_API_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/test/API
  -- HDF5_TEST_API_PAR_SRC_DIR = D:/work/hdf5_work/hdf5-1.14.3/testpar/API
  -- HDF5_PACKAGE_STRING = HDF5 1.14.3
  -- HDF5_PACKAGE_TARNAME = hdf5
  -- HDF5_PACKAGE_URL = http://www.hdfgroup.org
  -- HDF5_PACKAGE_BUGREPORT = help@hdfgroup.org
  -- HDF5_VERSION_STRING = 1.14.3
  -- HDF5_VERSION_MAJOR = 1.14
  -- HDF5_VERSION_MINOR = 3
  -- HDF_RESOURCES_DIR = D:/work/hdf5_work/hdf5-1.14.3/config/cmake
  -- HDF5_PACKAGE_NAME = HDF5  
  
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
  



