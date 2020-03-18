[Toc]
[English](./README.md) [简体中文 ](./README_zh_CN.md)
# OneFLOW

| Branch  | Linux Build | Windows Build |
|---      |---    |---    |
| Master  | [![Build Status](https://travis-ci.org/eric2003/OneFLOW.svg?branch=master)](https://travis-ci.org/eric2003/OneFLOW) | [![Build status](https://ci.appveyor.com/api/projects/status/o7fc231lp9jxlsib/branch/master?svg=true)](https://ci.appveyor.com/project/eric2003/OneFLOW/branch/master)  |
-----------------------------------------------------------
The Open-Source CFD Code
-----------------------------------------------------------
LargeScale Multiphysics Scientific Simulation Environment
-----------------------------------------------------------

## Build
Firstly, make sure that c++ compiler has been installed. For Windows platform, the [Visual Studio IDE](https://visualstudio.microsoft.com/ "Visual Studio IDE") is recommended. Then download and install [Git](https://git-scm.com/ "Git") and [Cmake](https://cmake.org/download/ "cmake") on your system.

### Dependencies

* [CGNS](https://github.com/CGNS/CGNS "CGNS"): a general, portable, and extensible standard for the storage and retrieval of CFD analysis data.
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/ "hdf5"): a set of file formats (HDF4, HDF5) designed to store and organize large amounts of data, needed by CGNS.
* [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download "Metis"): a set of serial programs for partitioning graphs, partitioning finite element meshes.
* [MPI](https://computing.llnl.gov/tutorials/mpi/ "MPI"): Massage Passing Interface.A standardized and portable message-passing standard for parallel computing. Some well-known open source implementations are [MS-MPI](https://github.com/Microsoft/Microsoft-MPI "MS-MPI")(recommended in Windows platform), [MPICH](https://github.com/pmodels/mpich "MPICH") and [OpenMPI](https://github.com/open-mpi/ompi "OpenMPI").

### Windows

1. Download source code from github:
```
git clone --recursive https://github.com/eric2003/OneFLOW
```
The above operation will download the source code together with prebuilt thirdparty libraries. You can also build them by yourself.

2. Use cmake to configure and generate project files.
   
3. Compile and generate executable file.
   
### Linux

1. Download source code from github:
```
git clone https://github.com/eric2003/OneFLOW
```
The above operation will download the source code only(prebuilt thirdparty libraries are not provided)

2. Use cmake to configure and generate appropriate project files, for example: Makefile.
   
3. Compile and generate executable file.
   
## OneFLOW DEVELOPERS
-----------------------------------------------------------
OneFLOW is being developed by individuals and organized teams all around the world.

The current OneFLOW release has been coordinated by the OneFLOW International Developers Society with selected contributions from the open-source community.

The main research teams contributing to the current release are:

Dr He Xin, at Transition technology

Email:<fantasy_2003_@hotmail.com>

If you have any problem in building or running the code, please do not hesitate to contact.

Copyright 2017-2020, He Xin, and the OneFLOW contributors.
