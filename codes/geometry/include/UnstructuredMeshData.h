/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#pragma once
#include "HXDefine.h"
#include <vector>
#include <string>
#include <map>
#include <set>

BeginNameSpace( ONEFLOW )

// Shared staging model for the supported first-order unstructured conversion
// path. Element codes currently use VTK identifiers and boundary types carry
// explicit CGNS mappings; this is not the complete canonical Mesh Core model.
const int MAX_VTK_TYPE = 100;
class VTK_TYPE
{
public:
    VTK_TYPE() = default;
    ~VTK_TYPE() = default;
public:
    static int VERTEX;
    static int LINE;
    static int TRIANGLE;
    static int QUADRILATERAL;
    static int TETRAHEDRON;
    static int HEXAHEDRON;
    static int PRISM;
    static int PYRAMID;
};

class VTK_CgnsMap
{
public:
    VTK_CgnsMap();
    ~VTK_CgnsMap();
public:
    std::map<int, int> vtk2Cgns;
    void Init();
};

extern VTK_CgnsMap vtk_CgnsMap;

class Marker
{
public:
    Marker() = default;
    ~Marker() = default;
public:
    std::string name;
    std::string bcName;
    int cgns_bcType{0};
    LinkField elems;
    IntField eTypes;
    int nElem{0};
};

class SecMarker
{
public:
    SecMarker() = default;
    ~SecMarker() = default;
public:
    int vtk_type{0};
    int cgns_type{0};
    std::string name;
    int nElem{0};
    LinkField elems;
};

class SecMarkerManager
{
public:
    SecMarkerManager() = default;
    ~SecMarkerManager() = default;
public:
    int nType{0};
    HXVector<SecMarker> data;
    void Alloc(int nType);
    int CalcTotalElem() const;
};

class MarkerManager
{
public:
    MarkerManager() = default;
    ~MarkerManager() = default;
public:
    int nMarker{0};
    HXVector<Marker> markerList;
    IntField types;
    LinkField l2g;
    void CreateMarkerList(int nMarker);
    void CalcSecMarker(SecMarkerManager &secMarkerManager);
};

class UnstructuredMeshData;

class VolumeSecManager
{
public:
    VolumeSecManager() = default;
    ~VolumeSecManager() = default;
public:
    int nSection{0};
    IntField types;
    LinkField l2g;
    void CalcVolSec(const UnstructuredMeshData &mesh, SecMarkerManager &secMarkerManager);
};

// Source readers populate this value model; CGNS-specific serialization is
// implemented by UnstructuredCgnsZoneWriter, not by a source-format reader.
class UnstructuredMeshData
{
public:
    UnstructuredMeshData();
    ~UnstructuredMeshData() = default;
public:
    int ndim{0};
    int nPoin{0};
    int nElem{0};
    IntField elemVTKType;
    LinkField elems;
    RealField xN;
    RealField yN;
    RealField zN;
    MarkerManager mmark;
    VolumeSecManager volSec;
};

EndNameSpace
