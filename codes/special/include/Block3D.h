/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "CalcCoor.h"
#include "BlkMesh.h"

BeginNameSpace( ONEFLOW )

class Block3D;

class MDomain;
class Face2D;
class Grid;
class StrGrid;

class Block3D : public BlkBasic
{
public:
    Block3D();
    ~Block3D();
public:
    RealField3D x3d, y3d, z3d;
    HXVector< MDomain * > mDomainList;
public:
    void Alloc();
    int GetNSubDomain();
    void ConstructTopo();
    void SetInterfaceBc();
    void GetCornerPoint( int & pt, int id1, int id2, int id3 );
    void CalcBlkDim();
    void CreateFaceList();
public:
    void CreateBlockMesh();
    void GenerateBlockMesh();
    void FillStrGrid( Grid * gridIn, int iZone );
};



EndNameSpace