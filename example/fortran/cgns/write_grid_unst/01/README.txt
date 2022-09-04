PS D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build> cmake ..
-- Building for: Visual Studio 17 2022
-- Selecting Windows SDK version 10.0.19041.0 to target Windows 10.0.22000.
-- The C compiler identification is MSVC 19.33.31629.0
-- The CXX compiler identification is MSVC 19.33.31629.0
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.33.31629/bin/Hostx64/x64/cl.exe - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Tools/MSVC/14.33.31629/bin/Hostx64/x64/cl.exe - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- The Fortran compiler identification is Intel 2021.5.0.20211109
-- Detecting Fortran compiler ABI info
-- Detecting Fortran compiler ABI info - done
-- Determine Intel Fortran Compiler Implicit Link Path
-- Determine Intel Fortran Compiler Implicit Link Path - done
-- Check for working Fortran compiler: C:/Program Files (x86)/Intel/oneAPI/compiler/2022.0.3/windows/bin/intel64/ifort.exe - skipped
-- The CGNS_INCLUDE_DIRS is C:/Program Files (x86)/cgns/include
-- Configuring done
-- Generating done
-- Build files have been written to: D:/work/fortran_work/ModernFortran/codes/cgns/write_grid_unst/01/build
PS D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build> cmake --build .
Build started...
1>------ Build started: Project: testprj, Configuration: Debug x64 ------
1>Compiling with Intel? Fortran Compiler Classic 2021.5.0 [Intel(R) 64]...
1>main.f90
1>Compiling manifest to resources...
1>Microsoft (R) Windows (R) Resource Compiler Version 10.0.10011.16384
1>Copyright (C) Microsoft Corporation.  All rights reserved.
1>Linking...
1>Microsoft (R) Incremental Linker Version 14.33.31629.0
1>Copyright (C) Microsoft Corporation.  All rights reserved.
1>/OUT:D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build\Debug\testprj.exe
1>/VERSION:0.0
1>/MANIFEST
1>/MANIFESTFILE:testprj.dir\Debug\testprj.exe.intermediate.manifest
1>"/MANIFESTUAC:level='asInvoker' uiAccess='false'"
1>/DEBUG
1>/PDB:D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build\Debug/testprj.pdb
1>/SUBSYSTEM:CONSOLE
1>/IMPLIB:D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build\Debug\testprj.lib
1>user32.lib
1>"C:\Program Files (x86)\cgns\lib\cgns.lib"
1>/machine:x64
1>/debug
1>/INCREMENTAL
1>testprj.dir\Debug\main.obj
1>testprj.dir\Debug\testprj.exe.embed.manifest.res
1>LINK : warning LNK4098: defaultlib 'MSVCRT' conflicts with use of other libs; use /NODEFAULTLIB:library
1>Embedding manifest...
1>Microsoft (R) Manifest Tool
1>Copyright (c) Microsoft Corporation.
1>All rights reserved.
1>
1>Build log written to  "file://D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build\testprj.dir\Debug\BuildLog.htm"
1>testprj - 0 error(s), 1 warning(s)
2>------ Build started: Project: ALL_BUILD, Configuration: Debug x64 ------
2>Building Custom Rule D:/work/fortran_work/ModernFortran/codes/cgns/write_grid_unst/01/CMakeLists.txt
========== Build: 2 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========
PS D:\work\fortran_work\ModernFortran\codes\cgns\write_grid_unst\01\build> .\Debug\testprj.exe
 Program write_grid_unst
 ...using 64-bit mode for particular integers
 created simple 3-D grid points
 Successfully wrote unstructured grid to file grid.cgns