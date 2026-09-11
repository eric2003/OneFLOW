CGNS
==================================

CGNS (Computational Fluid Dynamics General Notation System) is a universal file format and library for storing and exchanging computational fluid dynamics (CFD) data.

The CGNS format is scalable and portable, making it suitable for various CFD applications and tools. It provides a standardized way to describe and transmit CFD grids, flow fields, and other related data, enabling data exchange and interoperability between different CFD software.

The CGNS format contains a hierarchical structure and supports multiple grid types (such as structured, unstructured, hybrid grids, etc.) and cell types (such as triangles, tetrahedrons, hexahedrons, etc.). It also supports multiple bases and zones, allowing for separate processing and storage of different flow fields.

There are also open-source tools and libraries available for CGNS format, including CGNS API (C, Fortran, and Python), HDF5 library (for storing and reading CGNS files), ParaView (for visualizing CGNS data), and more.

Overall, CGNS format and library are important and useful tools in the field of CFD, simplifying and accelerating the storage, transmission, and processing of CFD data.


#. `CGNS and HDF5 compiling and linking related article links compilation <https://zhuanlan.zhihu.com/p/452874893/>`_
#. `CGNS official website <https://cgns.github.io/>`_
#. `CGNS standard document <https://cgns.github.io/CGNS_docs_current/>`_
#. `CGNS code repository <https://github.com/CGNS/CGNS/>`_
#. `CGNS4.4.0源码下载地址 <https://zhuanlan.zhihu.com/p/661758085/>`_
#. `CGNS从入门到精通系列链接整理 <https://zhuanlan.zhihu.com/p/560341685/>`_
#. `unst_example <https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_example/>`_


::

  cmake -DCMAKE_INSTALL_PREFIX=c:/dev/cgns/4.4.0 .. `
        -DCGNS_ENABLE_HDF5=ON `
        -DCGNS_BUILD_SHARED=ON `
        -DCGNS_USE_SHARED=ON `
        -DCGNS_BUILD_TESTING=ON `
        -DCGNS_ENABLE_TESTS=ON `
        -DCGNS_BUILD_CGNSTOOLS=ON `
        -DHDF5_DIR="c:/dev/HDF_Group/HDF5/1.14.2/cmake"

::

  cd d:\work\cgns_work\CGNS-4.4.0\build\

set oneAPI environment
::

  cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
  
Compile CGNS with HDF5
::

  cmake -DCMAKE_INSTALL_PREFIX=c:/dev/cgns/4.4.0 .. `
        -DCMAKE_BUILD_TYPE=Release `
        -DCGNS_ENABLE_HDF5=ON `
        -DCGNS_ENABLE_64BIT=ON `
        -DCGNS_ENABLE_FORTRAN=ON `
        -DCGNS_ENABLE_TESTS=ON `
        -DCGNS_BUILD_SHARED=ON `
        -DCGNS_USE_SHARED=ON `
        -DCGNS_BUILD_TESTING=ON `
        -DCGNS_BUILD_CGNSTOOLS=ON `
        -DHDF5_DIR="c:/dev/HDF_Group/HDF5/1.14.2/cmake"
  cmake --build . --parallel 4 --config Release
  cmake --install . --config Release --prefix c:/dev/cgns/4.4.0

Compile CGNS without HDF5
::

  cmake -DCMAKE_INSTALL_PREFIX=c:/dev/cgns/4.4.0 .. `
        -DCMAKE_BUILD_TYPE=Release `
        -DCGNS_ENABLE_HDF5=OFF `
        -DCGNS_ENABLE_64BIT=ON `
        -DCGNS_ENABLE_FORTRAN=ON `
        -DCGNS_ENABLE_TESTS=OFF `
        -DCGNS_BUILD_SHARED=OFF `
        -DCGNS_USE_SHARED=OFF

Compile CGNS without HDF5 BUILD_SHARED
::

  cmake -DCMAKE_INSTALL_PREFIX=c:/dev/cgns/4.4.0 .. `
        -DCMAKE_BUILD_TYPE=Release `
        -DCGNS_ENABLE_HDF5=OFF `
        -DCGNS_ENABLE_64BIT=ON `
        -DCGNS_ENABLE_FORTRAN=ON `
        -DCGNS_ENABLE_TESTS=OFF `
        -DCGNS_BUILD_SHARED=ON `
        -DCGNS_USE_SHARED=OFF
		
Modify CGNS CMake for mixed Fortran C
::

  #add_library(cgns_shared SHARED ${cgns_FILES} $<$<BOOL:${CGNS_ENABLE_FORTRAN}>:$<TARGET_OBJECTS:cgns_f2c>>)
  add_library(cgns_shared SHARED ${cgns_FILES} )
  target_link_directories( cgns_shared PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} )
  target_link_libraries(cgns_shared 
      PRIVATE $<$<BOOL:${CGNS_ENABLE_FORTRAN}>:$<TARGET_OBJECTS:cgns_f2c>>
  )


Start cmake gui from the command line
::

  & "C:\Program Files\CMake\bin\cmake-gui.exe"

powershell add environmental path
::

  $env:path += ";C:/dev/cgns/4.4.0/bin/"
  $env:path += ";C:/dev/HDF_Group/HDF5/1.14.2/bin/"
  
cgnsview
::

  modify
  proc file_stats {} {
    global ProgData
    set ProgData(file,size) "[file size $ProgData(file,name)] bytes"  
  to
  proc file_stats {} {
    global ProgData
    set ProgData(file,size) "[file size $ProgData(file,name)] bytes"
    if { [ catch {CGIOopen $ProgData(file)} msg] } {
        return $msg
    }
  
#. `unst_tetra  <https://cgns.github.io/CGNS_docs_current/sids/conv.html#unst_tetra/>`_
 
  
For all element types except MIXED, ElementConnectivity contains the list of nodes for each element. If the elements are sorted, then it must first list the connectivity of the boundary elements, then that of the interior elements.
::

  ElementConnectivity = Node11, Node21, ... NodeN1,
                        Node12, Node22, ... NodeN2,
                        ...
                        Node1M, Node2M, ... NodeNM  

where M is the total number of elements (i.e., ElementSize), and N is the number of nodes per element.
ElementDataSize indicates the total size (number of integers) of the array ElementConnectivity. For all element types except MIXED, NGON_n, and NFACE_n, ElementDataSize is given by:
::

  ElementDataSize = ElementSize * NPE[ElementType]
  
where NPE[ElementType] is a function returning the number of nodes for the given ElementType. For example, NPE[HEXA_8]=8.
When the section ElementType is MIXED, the data array ElementConnectivity contains one extra integer per element, to hold each individual element type:
::

  ElementConnectivity = Etype1, Node11, Node21, ... NodeN1,
                        Etype2, Node12, Node22, ... NodeN2,
                        ...
                        EtypeM, Node1M, Node2M, ... NodeNM						
						
where again M is the total number of elements, and Ni is the number of nodes in element i. The data array ElementStartOffset contains the starting positions of each element in the ElementConnectivity data array and its last value corresponds to the ElementConnectivity total size:
::

  ElementStartOffset  = 0, NPE[Etype1] + 1, ... ElementStartOffset[n-1] + NPE[Etypen] + 1,
                        ..., ElementStartOffset[M-1] + NPE[EtypeM] + 1 = ElementDataSize
						
In the case of MIXED element section, ElementDataSize is given by:
::

  ElementDataSize = ∑(NPE[ElementTypen] + 1)
  
where the summation is over n, and n represents a specific element type.
Arbitrary polyhedral elements may be defined using the NGON_n and NFACE_n element types. The NGON_n element type is used to specify all the faces in the grid, and the NFACE_n element type is then used to define the polyhedral elements as a collection of these faces. Except for boundary faces, each face of a polyhedral element must be shared by another polyhedral element.

For example, for NGON_n, the data array ElementConnectivity contains a list of nodes making up each face in the grid while ElementStartOffset provides the starting position of each face in the ElementConnectivity array:
::

  ElementConnectivity = Node11, Node21, ... NodeN1,
                        Node12, Node22, ... NodeN2,
                        ...
                        Node1M, Node2M, ... NodeNM

  ElementStartOffset  = 0, Nnodes1, Nnodes1 + Nnodes2, ...
                        ..., ElementStartOffset[i-1] + Nnodesi,
                        ..., ElementStartOffset[M-1] + NnodesM = ElementDataSize
						
where here M is the total number of faces, and Nnodesi is the number of nodes in face i. The ElementDataSize is the total number of nodes defining all the faces. Note that the number of nodes in face i is given by:
::

  Nnodesi = ElementStartOffset[i+1] - ElementStartOffset[i]
  
Then for NFACE_n, ElementConnectivity contains the list of face elements making up each polyhedral element, while ElementStartOffset provides the starting position of each polyhedral element in the ElementConnectivity array:
::

  ElementConnectivity = Face11, Face21, ... FaceN1,
                        Face12, Face22, ... FaceN2,
                        ...
                        Face1M, Face2M, ... FaceNM

  ElementStartOffset  = 0, Nfaces1, Nfaces1 + Nfaces2, ...
                        ..., ElementStartOffset[i-1] + Nfacesi,
                        ..., ElementStartOffset[M-1] + NfacesM = ElementDataSize						



