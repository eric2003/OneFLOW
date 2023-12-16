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



