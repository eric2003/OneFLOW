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
#include "CalcCoor.h"
#include "BlkMesh.h"

BeginNameSpace( ONEFLOW )

class Block3D;

class MDomain;
class Face2D;
class Grid;
class StrGrid;


class MLine;

class Block2D : public BlkBasic
{
public:
    Block2D();
    ~Block2D();
public:
    RealField2D x2d, y2d, z2d;
    HXVector< MLine * > mLineList;
    HXVector< MDomain * > mDomainList;
public:
    void Alloc();
    void CreateBlockMesh2D();
    void GenerateBlockMesh2D();
    int GetNSubDomain();
    void ConstructTopo();
    void SetInterfaceBc();
    void GetCornerPoint( int & pt, int id1, int id2 );
    void CalcBlkDim();
    void CreateFaceList();
    void FillStrGrid( Grid * gridIn, int iZone );
    void DumpBlockMesh2D( fstream &file );
};

EndNameSpace