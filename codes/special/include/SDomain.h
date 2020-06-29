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
#include "Mid.h"
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include <set>
#include <map>
#include <fstream>
using namespace std;
BeginNameSpace( ONEFLOW )

class Block3D;
class BlkMesh;
class SDomain;
class MDomain;
class CalcCoor;
class Face2D;
class MLine;

class SDomain : public DomData
{
public:
    SDomain();
    ~SDomain();
public:
    int domain_id;
    HXVector< MLine * > mLineList;
    CoorMap * localCoorMap;
    CoorMap * coorMap;
    RealField2D x2d, y2d, z2d;
public:
    void Alloc();
    void SetDomain( int fid, IntField & lineList, IntField & posList );
    void ConstructSDomainCtrlPoint();
    void GetCommonPoint( MLine * mLine1, MLine * mLine2, int & pt );
    bool CalcSingleDomainCoor();
    void SetRemainingCtrlPoint( IntField & idxList );
public:
    void ConstructLineToDomainMap();
    void ConstructLineToDomainMap( map< int, IntSet > & lineToDomainMap );
    void ConstructPointToDomainMap();
    void ConstructPointToDomainMap( map< int, IntSet > & pointToDomainMap );
    void ConstructPointToPointMap();
    void ConstructPointToPointMap( map< int, IntSet > & pointToPointMap );
    void ConstructDomainTopo();
    void GetPointIdLink( IntField & lineList, LinkField & pointIdLink );
public:
    void Add( IntField &iList, IntField &jList, IntField &kList, int i, int j, int k );
    void ConstructLocalTopoAsBlk2D();
    void SetBlkBcMesh( Block3D * blk3d );
    void SetDomainBcMesh();
    void GenerateSDomainMesh();
    void GenerateSDomainMesh( fstream & file );
};


EndNameSpace