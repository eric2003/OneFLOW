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
#include "HXCgns.h"
#include "Mid.h"
#include <vector>
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )
class PointFactory;
class ElemFeature;
class FaceTopo;

class FaceSolver
{
public:
    FaceSolver();
    ~FaceSolver();
public:
    set< Mid<int> > * refFaces;
    IntField * faceBcKey;
    IntField * faceBcType;
    LinkField * childFid;
public:
    FaceTopo * faceTopo;
public:
    int FindFace( Mid<int> & face );
    bool CheckBcFace( IntSet & bcVertex, IntField & nodeId );
    void ScanElementFace( CgIntField & eNodeId, int eType, int eId );
    void ScanBcFace( IntSet & bcVertex, int bcType, int bcNameId );
    void ScanBcFaceDetail( IntSet & bcVertex, int bcType, int bcNameId );
    void ScanInterfaceBc();
    int GetNSimpleFace();
};

EndNameSpace