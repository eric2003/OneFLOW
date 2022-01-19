/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class Block3D;

class MDomain;
class Face2D;
class Grid;
class StrGrid;

class BlkBasic
{
public:
    BlkBasic();
    virtual ~BlkBasic();
public:
    int blk_id;
    int ni, nj, nk;
    LinkField localpt;
    IntField controlpoints;
    HXVector< Face2D * > facelist;
    IntField dimList;
    CoorMap coorMap;
public:
    void AddLocalPt( int p1, int p2 );
    void AddLocalPt( int p1, int p2, int p3, int p4 );
    void DumpInp( std::fstream & file );
    virtual int GetNSubDomain() { return 1;  };
public:
    void Add( IntField &iList, IntField &jList, IntField &kList, int i, int j, int k );
};

void SetGridXYZ3D( StrGrid * grid, RealField3D & x3d, RealField3D & y3d, RealField3D & z3d );
void SetGridXYZ2D( StrGrid * grid, RealField2D & x2d, RealField2D & y2d, RealField2D & z2d );


EndNameSpace
