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
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class InterFace;

class BcInfo
{
public:
    BcInfo();
    ~BcInfo();
public:
    LinkField bcFace;
    LinkField bcRegion;
    LinkField bcdtkey;
    IntField bcType;
public:
    UInt GetNBcRegion() { return bcType.size(); }
};

class BcRecord
{
public:
    BcRecord();
    ~BcRecord();
public:
    IntField bcType;
    IntField bcdtkey;
    IntField bcRegion;
    BcInfo * bcInfo;
public:
    void Init( UInt nBFace );
    int GetNBFace();
    int CalcNIFace();
    int CalcNumWallFace();
    void CreateI2B( InterFace * interFace );
public:
    void CreateBcRegion();
};

class IFaceLink;

class BcManager
{
public:
    BcManager();
    ~BcManager();
public:
    bool deleteBoundaryCondition;
    BcRecord * bcRecord;
    BcRecord * bcRecordNew;
    IntField l2gNew;

    IntField bcKeyVector;

    LinkField childFaceId;
    IntField numberOfChildFacesOfInterface;
    IntField newBoundaryLeftCellIndex;

    IntField bcFlag;
public:
    void PreProcess();
    bool ExistInterface();
    void Update();
    void CalcBcType( IntField & bcTypeList );
};

class BasicRegion
{
public:
    int zid;
    int dir;         //边界面方向:0, 1, 2对应于i, j, k
    int outerNormal; //左右边界-1, 1对应于左右边界
    int start[ 3 ], end[ 3 ];
    int lr[ 3 ];     //左右边界-1, 1对应于左右边界
public:
    void SetRegion( int ist, int ied, int jst, int jed );
    void SetRegion( int ist, int ied, int jst, int jed, int kst, int ked );
};

class BasicRegion;
class TestRegion
{
public:
    TestRegion();
    ~TestRegion();
    int p1[ 3 ], p2[ 3 ];
    int a[ 3 ];
    int sign[ 3 ];
public:
    void Run( BasicRegion * r, int dimension );
};

class BcRegion;
class TestRegionM
{
public:
    TestRegionM();
    ~TestRegionM();
public:
    int dimension;
    TestRegion s, t;
    int itransform[ 3 ];
public:
    void Run( BcRegion * bcRegion, int dimension );
};

class BcRegion
{
public:
    BcRegion( int zid, int rid );
    ~BcRegion();
public:
    int rid;                         //region id
    int bcType;                      //boundary type
    string regionName;               //boundary name
public:
    BasicRegion * s;
    BasicRegion * t;
public:
    void GetNormalizeIJKRegion( int & ist, int & ied, int & jst, int & jed, int & kst, int & ked );
    int CalcRegionCells();
};

class BcRegionGroup
{
public:
    BcRegionGroup();
    ~BcRegionGroup();
public:
    int zoneIndex;
    int nBFace, nIFace;
    HXVector< BcRegion * > * regions;
    void Create( int nBcRegions );
    void SetBcRegion( int ir, BcRegion * bcRegion );
    BcRegion * GetBcRegion( int ir );
};

EndNameSpace