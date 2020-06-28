/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
using namespace std;

BeginNameSpace( ONEFLOW )
class FaceTopo;
class NodeMesh;
class CellMesh;

class FaceMesh
{
public:
    FaceMesh();
    ~FaceMesh();
public:
    FaceTopo * faceTopo;
    RealField xfc, yfc, zfc;
    RealField xfn, yfn, zfn;
    RealField area;
    RealField vfx, vfy, vfz;
    RealField vfn;
public:
    UInt GetNFace();
    UInt CalcTotalFaceNodes();
    UInt GetNBFace();
    void SetNBFace( UInt nBFace );
    void CalcFaceNormal1D( NodeMesh * nodeMesh, CellMesh * cellMesh );
    void CalcFaceCenter1D( NodeMesh * nodeMesh );
    void CalcFaceNormal2D( NodeMesh * nodeMesh );
    void CalcFaceCenter2D( NodeMesh * nodeMesh );
    void CalcFaceNormal3D( NodeMesh * nodeMesh );
    void CalcFaceCenter3D( NodeMesh * nodeMesh );

    void AllocateMetrics();
};


EndNameSpace