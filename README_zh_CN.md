[Toc]
[English](./README.md) [简体中文 ](./README_zh_CN.md)
# OneFLOW
-----------------------------------------------------------
The Open-Source CFD Code
-----------------------------------------------------------
LargeScale Multiphysics Scientific Simulation Environment
-----------------------------------------------------------

## 编译
首先,确认系统中已经安装了C++编译器. 对于Windows平台, 建议使用[Visual Studio IDE](https://visualstudio.microsoft.com/ "Visual Studio IDE")。然后下载并安装[Git](https://git-scm.com/ "Git")和[Cmake](https://cmake.org/download/ "cmake")。

### 第三方依赖库

* [CGNS](https://github.com/CGNS/CGNS "CGNS"): 用于处理CFD网格和数据的一项标准。
* [HDF5](https://www.hdfgroup.org/downloads/hdf5/ "hdf5"): 用于处理海量数据的一种文件存储标准, 编译CGNS时需要使用HDF5。
* [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/download "Metis"): 用于对图结构和网格进行分区的一个程序库。
* [MPI](https://computing.llnl.gov/tutorials/mpi/ "MPI"): 消息传递接口，一个可移植，高性能的并行计算标准. 开源的实现有 [MS-MPI](https://github.com/Microsoft/Microsoft-MPI "MS-MPI")(适用于Windows平台), [MPICH](https://github.com/pmodels/mpich "MPICH")和[OpenMPI](https://github.com/open-mpi/ompi "OpenMPI")。

库的编译：
    CGNS，HDF5和Metis均使用Cmake进行编译管理，MS-MPI可直接下载编译好的库。

### Windows

1. 从github下载源代码:
```
git clone --recursive https://github.com/eric2003/OneFLOW
```
上述操作将会下载源代码及编译好的第三方依赖库，也可以自行编译这些库。

2. 使用Cmake设置编译选项，配置并生成相应的项目文件。如果采用Visual Studio进行编译，并且下载了第三方库，则Cmake会自动配置库的路径。否则需要在Cmake的配置界面修改相应的默认设置。
   
3. 编译并生成可执行文件。
   
### Linux

1. 下载或者编译第三方依赖库
   debian/ubuntu系统：
   ```
   sudo apt install libhdf5-dev libcgns-dev libmetis-dev libmpich-dev
   ```

   以ubuntu16.04为例，使用上述命令安装第三方库后，各个库的相应信息为：

    cgns  

    头文件路径：/usr/include

    库文件路径：/usr/lib  

    库文件名称：cgns

    hdf5 

    头文件路径：/usr/include/hdf5

    库文件路径：/usr/lib/x86_64-linux-gnu/hdf5/serial

    库文件名称：hdf5

    metis

    头文件路径：/usr/include

    库文件路径：/usr/lib/x86_64-linux-gnu

    库文件名称：metis

    mpich 

    头文件路径：/usr/include/mpich

    库文件路径：/usr/lib

    库文件名称：mpi

2. 从github下载源代码:
```
git clone https://github.com/eric2003/OneFLOW
```
上述操作将会下载源代码（不包括编译好的第三方依赖库）

3. 使用Cmake设置编译选项，配置并生成相应的项目文件, 如Linux系统常用的Makefile文件。
   
1. 编译并生成可执行文件。
   
## OneFLOW开发者
-----------------------------------------------------------
OneFLOW由分散的团队和个人共同开发。

The current OneFLOW release has been coordinated by the OneFLOW International Developers Society with selected contributions from the open-source community.

当前代码的主要贡献团队:

赫新 博士, 转捩点科技

电子邮箱：<fantasy_2003_@hotmail.com>

如果在编译和运行代码中遇到问题，可随时通过邮件联系。

Copyright 2017-2020, He Xin, and the OneFLOW contributors.
