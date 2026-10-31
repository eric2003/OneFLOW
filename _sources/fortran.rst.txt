Fortran
==================================

Fortran, as derived from Formula Translating System, is a general-purpose, imperative programming language. It is used for numeric and scientific computing.

#. `Fortran Programming Language <https://fortran-lang.org/>`_
#. `Fortran Tutorial <https://www.tutorialspoint.com/fortran/index.htm>`_
#. `Getting started with gfortran (PDF) <https://gcc.gnu.org/onlinedocs/gfortran.pdf>`_
#. `gfortran(1) — Linux manual page <https://man7.org/linux/man-pages/man1/gfortran.1.html>`_
#. `Using GNU Fortran <https://gcc.gnu.org/wiki/GFortranGettingStarted/>`_
#. `INTRODUCTION TO FORTRAN <https://ourcodingclub.github.io/tutorials/fortran-intro/>`_
#. `GNU Fortran compiler (gfortran) <http://magnin.plil.net/spip.php?article97>`_
#. `LINK : fatal error LNK1104: cannot open file 'ifmodintr.lib' <https://community.intel.com/t5/Intel-Fortran-Compiler/LINK-fatal-error-LNK1104-cannot-open-file-ifmodintr-lib/td-p/933371/>`_
#. `Configuring Visual Studio for Mixed-Language Applications <https://www.intel.com/content/www/us/en/developer/articles/technical/configuring-visual-studio-for-mixed-language-applications.html>`_
#. `Fortran for C/C++ developers made easier with CMake <https://www.kitware.com/fortran-for-cc-developers-made-easier-with-cmake/>`_
#. `Using C/C++ and Fortran together <http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html>`_
#. `Calling Fortran From C - Part 1: Hello, World <https://www.youtube.com/watch?v=urcy6-kXZDw/>`_
#. `Calling Fortran From C - Part 2: A Simple Function <https://www.youtube.com/watch?v=LmVUTWQDdC4/>`_
#. `Calling Fortran From C - Part 3: Arrays and Matrices <https://www.youtube.com/watch?v=_h8eZ7vI_uw/>`_
#. `Fortran与C/C++混合编程示例 <https://www.cnblogs.com/snake553/p/6962386.html>`_
#. `C与Fortran混合编程 <https://zhejianggaoxiao.github.io/2017/03/15/C%E4%B8%8EFortran%E6%B7%B7%E5%90%88%E7%BC%96%E7%A8%8B/>`_
#. `Visual Studio 中 Fortran 和 C 的混合编程 <https://zhuanlan.zhihu.com/p/455728972/>`_
#. `Mix C Fortran <https://wiki.ubuntu.org.cn/Mix_C_Fortran/>`_
#. `C-Fortran 接口 <https://docs.oracle.com/cd/E19205-01/820-1204/6nct259sc/index.html>`_
#. `Fortran-C-CPP混合编程-1 <https://www.jianshu.com/p/aefef17cef69/>`_
#. `混合 Fortran 和 C++ <https://www.ibm.com/docs/zh/openxl-fortran-aix/17.1.1?topic=calls-mixing-fortran-c/>`_
#. `Windows系统下Fortran编程 <https://www.bilibili.com/video/BV1XD4y1S7jz/>`_
#. `VS2019中C++与Fortran的混合编程 <https://blog.csdn.net/weixin_41124748/article/details/120037882/>`_
#. `Fortran C/C++ interoperability <https://github.com/aerosayan/fortran-c-interop/>`_
#. `A Modern Fortran Scientific Programming Ecosystem <https://degenerateconic.com/a-modern-fortran-scientific-programming-ecosystem.html>`_
#. `Intel® C++ & Fortran Compiler <https://indico.cern.ch/event/403113/contributions/1847270/attachments/1123921/1603880/03_Intel_C__Fortran_Compiler.pdf>`_
#. `Fortran Programming Tutorials (Revised) <https://www.youtube.com/watch?v=ProXdx6xJb8/>`_
#. `Compiling C/C++/Fortran code <https://wiki.usask.ca/pages/viewpage.action?pageId=1955337523>`_
#. `Fortran with CMake (Simple Tutorial) <https://www.youtube.com/watch?v=Tl3Ph-4dMTI/>`_
#. `Modern Fortran: Concurrency and Parallelism <https://www.youtube.com/watch?v=tdjo1OI-30g/>`_
#. `Parallel programming without MPI - Using coarrays in Fortran <https://www.youtube.com/watch?v=tGSoCvTLfkw/>`_
#. `Fortran 90 Module Dependencies <https://lagrange.mechse.illinois.edu/f90_mod_deps/>`_
#. `Modern Fortran in Science and Technology <https://modern-fortran-in-science-and-technology.readthedocs.io/en/latest/index.html>`_
#. `Writing Makefiles for Modern Fortran <https://aoterodelaroza.github.io/devnotes/modern-fortran-makefiles/>`_
#. `Mixing C++ and Fortran <https://enccs.github.io/cmake-workshop/cxx-fortran/>`_
#. `Getting started with Fortran <https://riptutorial.com/fortran/>`_
#. `Visual Studio Code C/C++/Fortran with Multiple Source Files <https://iraspa.org/blog/visual-studio-code-c-cpp-fortran-with-multiple-source-files/>`_
#. `Using gfortran with external libraries and module files <https://kacv.net/brad/engr325/gfortrancompile.pdf>`_
#. `GfortranApps <https://gcc.gnu.org/wiki/GfortranApps>`_
#. `Setting Up Windows For Fortran Development <https://www.youtube.com/watch?v=kzloL99wtN0/>`_
#. `Awesome Fortran <https://github.com/rabbiabram/awesome-fortran/>`_



Compiling the source code
-------------------------------
But the basics are simple enough. Take the gfortran compiler, part of the GNU compiler collection. To compile a simple program as the one above, that consists of one source file, you run the following command, assuming the source code is stored in the file “hello.f90”:

::

  gfortran -c hello.f90
  
This results in a file “hello.o” (as the gfortran compiler uses “.o” as the extension for the object files).

The option “-c” means: only compile the source files. If you were to leave it out, then the default action of the compiler is to compile the source file and start the linker to build the actual executable program. The command:  

::

  gfortran hello.f90

::

  gfortran hello.f90 -o main
  
Building Dynamic-Link Libraries
----------------------------------
#. `Building Dynamic-Link Libraries <https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-reference-build-windows-applications/15-0/building-dynamic-link-libraries.html>`_

There is one more thing to be aware of: On Windows you must explicitly specify that a procedure is to be exported, i.e. is visible in the dynamic library. There are several ways — depending on the compiler you use — to achieve this. One method is via a so-called compiler directive:
::

  subroutine myroutine( ... )
  !GCC$ ATTRIBUTES DLLEXPORT:: myroutine
  
Or, with the Intel Fortran compiler:  
::
  
  subroutine myroutine( ... )
  !DEC$ ATTRIBUTES DLLEXPORT:: myroutine

Sample program: C calling Fortran

#. `Sample program: C calling Fortran <https://www.ibm.com/docs/en/xl-c-aix/13.1.0?topic=fortran-sample-program-calling/>`_
#. `Compiling a mixed C-Fortran program (main program is Fortran) <https://gcc.gnu.org/wiki/GFortranGettingStarted>`_
#. `C-Fortran Interface- <https://w.astro.berkeley.edu/~wright/f2c.html>`_
#. `FortranCon2020 [SP]: Shroud: generate Fortran wrappers for C and C++ libraries <https://www.youtube.com/watch?v=1mdI-M94vDc/>`_


::

  cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
  
Temporarily change the value of the environment variable PATH
::

  $env:path += ";d:\work\fortran_work\ModernFortran\codes\windows\shared-lib\01\"  
  
  
Building Shared Libraries
---------------------------------
#. `A.1.1 Building Shared Libraries <https://manuals.dianafea.com/d102/Analys/node622.html>`_
#. `Creating FORTRAN Libraries <https://www1.udel.edu/topics_css/software/special/language/fortran/fortran-docs/fortran-libraries1.html>`_

Fortran Formats
---------------------------
#. `Fortran Formats <https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html>`_

Mixed-Programming
------------------------------
#. `CALLING FORTRAN SUBROUTINES FROM FORTRAN, C, AND C++ <https://faculty.sites.iastate.edu/keinert/files/inline-files/calling.pdf>`_
#. `C-Fortran Interface <https://docs.oracle.com/cd/E19422-01/819-3685/11_cfort.html>`_
#. `新语法系列 之 interface 功能详解 <http://fcode.cn/guide-61-1.html>`_
#. `Calling 'C' from FORTRAN <https://community.intel.com/t5/Intel-Fortran-Compiler/Calling-C-from-FORTRAN/td-p/985459>`_
#. `Calling C++ from Fortran <https://fortran-lang.discourse.group/t/calling-c-from-fortran/3402>`_
#. `Jean Zay: Calling C functions from Fortran <http://www.idris.fr/eng/jean-zay/cpu/jean-zay-cpu-fortran_c-eng.html>`_
#. `Writing a Fortran-C interface <https://hackmd.io/@python-fortran-interface/rySynUcYd>`_
#. `Interoperation of Fortran with C <https://doku.lrz.de/files/10746213/10746217/1/1684600341047/Advanced_Fortran_Interop.pdf>`_
#. `ISO_C_BINDING <https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html>`_
#. `Iso_c_binding: Looking for practical example of how it helps with mangling <https://fortran-lang.discourse.group/t/iso-c-binding-looking-for-practical-example-of-how-it-helps-with-mangling/3393>`_
#. `Interoperability : Calling Fortran from C (i) <https://craftofcoding.wordpress.com/2018/01/29/interoperability-calling-fortran-from-c-i/>`_ 
#. `Mixing C++ and Fortran <https://enccs.github.io/cmake-workshop/cxx-fortran/>`_
#. `Philip Semanchuk| Python, C, C++, and Fortran Relationship Status: It’s Not That Complicated <https://www.youtube.com/watch?v=aUSokzzsEko/>`_
#. `Mixing C, C++, and Fortran <https://indico.ictp.it/event/a13229/session/2/contribution/11/material/0/0.pdf>`_
#. `Interoperable-Subroutines-and-Functions <https://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html>`_
#. `Interoperation of Fortran with C <https://doku.lrz.de/files/10746213/10746217/1/1684600341047/Advanced_Fortran_Interop.pdf>`_
#. `Iso_c_binding: Looking for practical example of how it helps with mangling <https://fortran-lang.discourse.group/t/iso-c-binding-looking-for-practical-example-of-how-it-helps-with-mangling/3393>`_
#. `Managing libraries (static and dynamic libraries) <https://fortran-lang.org/learn/building_programs/managing_libraries/>`_
#. `More on Calling C from FORTRAN <http://www.starlink.ac.uk/docs/sun209.htx/sun209se5.html>`_
#. `Mixed C/Fortran Programming <https://abhila.sh/writing/1.0/cfmix.html>`_
#. `C-Fortran Interface <https://docs.oracle.com/cd/E19422-01/819-3685/11_cfort.html>`_
#. `C structs in Fortran <https://riptutorial.com/fortran/example/7150/c-structs-in-fortran/>`_
#. `Using Fortran from C/C++ <https://github.com/wallytutor/medium-articles/tree/main/medium/2023-05-25-Calling-Fortran-From-C%2B%2B/interoperability-with-c>`_

fortran2018-examples
-------------------------
#. `fortran2018-examples <https://github.com/scivision/fortran2018-examples/>`_
#. `fortran-cpp-interface <https://github.com/scivision/fortran-cpp-interface/>`_
#. `fortran-submodule <https://github.com/scivision/fortran-submodule/>`_
#. `Fortran MPI Examples <https://github.com/scivision/fortran-mpi-examples/>`_
#. `Fortran Parallel Examples <https://github.com/scivision/fortran-parallel-examples/>`_
#. `Sparse Fortran libraries <https://github.com/scivision/sparse-fortran/>`_
#. `Object-oriented Fortran NetCDF4 interface <https://github.com/geospace-code/nc4fortran/>`_
#. `Geospace code <https://github.com/geospace-code/>`_


FORTRAN系列链接整理(FORTRAN series link)
--------------------------------------------
#. `FORTRAN系列链接整理 <https://zhuanlan.zhihu.com/p/662370470/>`_
#. `ubuntu22.04下查看gfortran版本号 <https://zhuanlan.zhihu.com/p/662369466/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-Hello, World示例代码 <https://zhuanlan.zhihu.com/p/662372522/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-Hello, World示例代码 version 1 <https://zhuanlan.zhihu.com/p/662383565/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-Hello, World示例代码 version 2 <https://zhuanlan.zhihu.com/p/662416747/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-addNumbers示例代码 <https://zhuanlan.zhihu.com/p/662423621/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-testingInt示例代码 <https://zhuanlan.zhihu.com/p/662434505/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-testingInt+integer+kind示例代码 <https://zhuanlan.zhihu.com/p/662436799/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-division+Real Type示例代码 <https://zhuanlan.zhihu.com/p/662447514/>`_
#. `ubuntu22.04+gfortran11.4.0编译运行Fortran-variableTesting示例代码 <https://zhuanlan.zhihu.com/p/662450844/>`_
#. `ubuntu+gfortran+创建并使用subroutine+static lib简单测试 <https://zhuanlan.zhihu.com/p/662955043/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+static lib简单测试 <https://zhuanlan.zhihu.com/p/664984017/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+static lib简单测试 <https://zhuanlan.zhihu.com/p/665024598/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665032940/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(module)+static lib简单测试 <https://zhuanlan.zhihu.com/p/665063802/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(module)+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665069959/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+files)+static lib简单测试 <https://zhuanlan.zhihu.com/p/665072636/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+files)+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665079819/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+files)+static lib v2简单测试 <https://zhuanlan.zhihu.com/p/665080732/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+files)+static lib简单测试 <https://zhuanlan.zhihu.com/p/665088630/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+files)+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665091815/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+files)+static lib v2简单测试 <https://zhuanlan.zhihu.com/p/665093599/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+subroutines+files)+static lib简单测试 <https://zhuanlan.zhihu.com/p/665097378/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+subroutines+files)+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665099647/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+subroutines+files)+static lib简单测试 <https://zhuanlan.zhihu.com/p/665102187/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+subroutines+files)+static lib v1 简单测试 <https://zhuanlan.zhihu.com/p/665103808/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+static lib v1简单测试 <https://zhuanlan.zhihu.com/p/665106327/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(single subroutine)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665123140/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(single subroutine)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665126451/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(single subroutine)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665136987/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(single module)+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665138079/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(single module)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665140915/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(single module)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665155494/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(single module)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665156269/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+files)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665157654/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+files)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665159392/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+files)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665176299/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+files)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665177428/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+subroutines+files)+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665179293/>`_
#. `windows11+oneAPI+icx+ifort+c调用Fortran+(multiple modules+subroutines+files)+shared lib v1 简单测试 <https://zhuanlan.zhihu.com/p/665181216/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+subroutines+files)+shared lib简单测试 <https://zhuanlan.zhihu.com/p/665182857/>`_
#. `ubuntu+gcc+gfortran+c调用Fortran+(multiple modules+subroutines+files)+shared lib v1简单测试 <https://zhuanlan.zhihu.com/p/665184053/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+say_hello+static lib 简单测试 <https://zhuanlan.zhihu.com/p/665236262/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+say_hello+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665252985/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+print_string+static lib 简单测试 <https://zhuanlan.zhihu.com/p/665289536/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+print_string+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665326639/>`_
#. `ubuntu+gcc+gfortran+Fortran调用c+print_string+static lib 简单测试 <https://zhuanlan.zhihu.com/p/665329906/>`_
#. `ubuntu+gcc+gfortran+Fortran调用c+print_string+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665333322/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+c_add_integer+static lib 简单测试 <https://zhuanlan.zhihu.com/p/665345253/>`_
#. `windows11+oneAPI+icx+ifort+Fortran调用c+c_add_integer+shared lib 简单测试 <https://zhuanlan.zhihu.com/p/665347386/>`_


main.c
::

  extern void show_N1();
  extern void show_N2();
  extern void show_N3();
  
  int main(int argc, char *argv[])
  {
      show_N1();
      show_N2();
      show_N3();
  
      return 0;
  }

onemod.f90
::

  module onemod  
  implicit none 
  
     integer, parameter :: N = 1024 
     
  contains      
     subroutine show_N() bind(C, name='show_N')
        print*, "N = ", N          
     end subroutine show_N 
     
  end module onemod 


Ubuntu+gcc+gfortran
::

  $ gfortran -c onemod.f90
  $ gfortran -c twomods.f90
  $ ar r mods.a onemod.o twomods.o
  ------------------------------------------
  check mods.a
  ------------------------------------------
  $ nm mods.a
  ------------------------------------------
  test(fortran)
  $ gfortran -c main.f90 -I../
  $ gfortran -o testprj main.o ../mods.a
  ------------------------------------------
  test(c)
  $ gcc -c main.c
  $ gcc -o testprj main.o ../mods.a -lgfortran
  
Windows11+oneAPI+icx+ifort
::

  $ cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
  
  $ ifort -c onemod.f90
  $ lib /OUT:onemod.lib onemod.obj
  ----------------------------------
  check onemod.lib
  ------------------------------------------
  $ dumpbin /symbols onemod.lib
  ------------------------------------------
  test(fortran)
  $ ifort -c main.f90
  $ ifort -o testprj main.obj ../onemod.lib
  -------------------------------------
  test(c)
  $ icx -c main.c
  $ icx -o testprj main.obj ../onemod.lib

Dumpbin
::

  PS D:\work\fortran_work\ModernFortran\codes\windows\c-call-fortran-lib\static\01> dumpbin /symbols sub.lib
  Microsoft (R) COFF/PE Dumper Version 14.37.32825.0
  Copyright (C) Microsoft Corporation.  All rights reserved.
  
  
  Dump of file sub.lib
  
  File Type: LIBRARY
  
  COFF SYMBOL TABLE
  000 00000000 SECT1  notype       Static       | .text
      Section length   50, #relocs    3, #linenums    0, checksum        0
  002 00000000 SECT1  notype ()    External     | sub_
  003 00000000 SECT2  notype       Static       | .rdata
      Section length   10, #relocs    0, #linenums    0, checksum        0
  005 00000008 SECT2  notype       Static       | __STRLITPACK_1
  006 00000000 UNDEF  notype ()    External     | for_write_seq_lis
  007 00000000 SECT3  notype       Static       | .xdata
      Section length    8, #relocs    0, #linenums    0, checksum        0
  009 00000000 SECT4  notype       Static       | .pdata
      Section length    C, #relocs    3, #linenums    0, checksum        0
  00B 00000000 SECT2  notype       Static       | __STRLITPACK_3.0.1
  00C 00000000 UNDEF  notype ()    External     | __ImageBase
  00D 00000000 SECT5  notype       Static       | .drectve
      Section length   B9, #relocs    0, #linenums    0, checksum        0
  
  String Table Size = 0x44 bytes
  
    Summary
  
            B9 .drectve
             C .pdata
            10 .rdata
            50 .text
             8 .xdata

Example1
--------------------------------------------

main.c
::

  extern void hello_print();
  int main( int argc, char**argv )
  {
  	hello_print();
  	return 0;
  }
  
sub.c
::

  #include <stdio.h>
  
  void hello_print()
  {
      printf("hello!\n");
  }

Windows11+oneAPI+icx
::

  cd d:\work\modern_cmake_work\ModernCMake\codes\cmake\add_library\test\07b\build\
  $ cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
  $ icx -c ../sub.c
  $ icx -c ../main.c 
  $ icx -o testprj main.obj sub.obj 
  
Windows11+oneAPI+icx static lib
::  

  cd d:\work\modern_cmake_work\ModernCMake\codes\cmake\add_library\test\07c\build\
  $ cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
  $ icx -c ../sub.c
  $ lib /OUT:sub.lib sub.obj
  $ icx -c ../main.c 
  $ icx -o testprj main.obj sub.lib 

  
Windows11+oneAPI+icx static lib+CMake
::  

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
  
  add_library(sub OBJECT sub.c)
  
  add_executable ( ${PROJECT_NAME}
      main.c
  )
  
  target_link_libraries( ${PROJECT_NAME} 
      PRIVATE 
          sub
  )
  
Windows11+oneAPI+icx static lib+CMake v1
::  

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
  
  add_library(sub OBJECT sub.c)
  
  add_executable ( ${PROJECT_NAME}
      main.c
  )
  
  target_link_libraries( ${PROJECT_NAME} 
      PRIVATE 
          $<TARGET_OBJECTS:sub>
  )
  
  get_target_property(prj_LINK_LIBRARIES ${PROJECT_NAME}  LINK_LIBRARIES)
  
  add_custom_target ( print 
      ${CMAKE_COMMAND} -E 
      echo 
      prj_LINK_LIBRARIES = ${prj_LINK_LIBRARIES} &
      echo 
      TARGET_OBJECTS:sub = $<TARGET_OBJECTS:sub>
  )
  
Windows11+oneAPI+icx static lib+CMake v2
::

  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
  
  add_library(sub OBJECT sub.c)
  
  add_executable ( ${PROJECT_NAME}
      main.c
      $<TARGET_OBJECTS:sub>
  )
  
  message ( STATUS "sub = ${sub}" )
  
  get_target_property(prj_LINK_LIBRARIES ${PROJECT_NAME}  LINK_LIBRARIES)
  
  add_custom_target ( print 
      ${CMAKE_COMMAND} -E 
      echo 
      prj_LINK_LIBRARIES = ${prj_LINK_LIBRARIES} &
      echo 
      TARGET_OBJECTS:sub = $<TARGET_OBJECTS:sub>
  )
  
::

  cmake ..
  cmake --build . --target print
  cmake --build . --config Debug --target print
  cmake --build . --config Release --target print  
  
Example2
--------------------------------------------
sub.c
::

  #include <stdio.h>
  
  void hello_print()
  {
      printf("hello!\n");
  } 

main.c
::

  extern void hello_print();
  int main( int argc, char**argv )
  {
  	hello_print();
  	return 0;
  }


fmain.f90
::

  program main
    implicit none
    
    interface
      subroutine fortran_hello_print() bind(C,name='hello_print')
      end subroutine fortran_hello_print
    end interface  
    
    call fortran_hello_print()
    
  end program main


Windows11+oneAPI+icx+ifort+CMake
::
  
  cmake_minimum_required ( VERSION 3.28 )
  
  project ( testprj )
  
  add_library(sub OBJECT sub.c)
  
  add_executable ( CPrj
      main.c
  )
  
  enable_language(Fortran)
  
  add_executable ( FPrj
      fmain.f90
  )
  
  target_link_libraries( CPrj
      PRIVATE 
          $<TARGET_OBJECTS:sub>
  )
  
  target_link_libraries( FPrj
      PRIVATE 
          $<TARGET_OBJECTS:sub>
  )
  
  get_target_property(cprj_LINK_LIBRARIES CPrj LINK_LIBRARIES)
  get_target_property(fprj_LINK_LIBRARIES FPrj LINK_LIBRARIES)
  
  add_custom_target ( print 
      ${CMAKE_COMMAND} -E 
      echo 
      cprj_LINK_LIBRARIES = ${cprj_LINK_LIBRARIES} &
      echo 
      fprj_LINK_LIBRARIES = ${fprj_LINK_LIBRARIES} &
      echo 
      TARGET_OBJECTS:sub = $<TARGET_OBJECTS:sub>
  )
  
 
  
::

  cmake ..
  cmake --build . --target print
  cmake --build . --config Debug --target print
  cmake --build . --config Release --target print
  
::

  Run Build Command(s): "C:/Program Files/Microsoft Visual Studio/2022/Community/Common7/IDE/devenv.com" IntelFortranImplicit.sln /build Debug /project ALL_BUILD  
  
::

  C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Visual Studio 2022\Visual Studio Tools
  C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Visual Studio 2022\Visual Studio Tools\Developer Command Prompt for VS 2022.lnk
  %comspec% /k "C:\Program Files\Microsoft Visual Studio\2022\Community\Common7\Tools\VsDevCmd.bat"
