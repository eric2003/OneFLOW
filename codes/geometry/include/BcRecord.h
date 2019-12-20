/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
    int ComputeNIFace();
    int CmpNumWallFace();
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
    void CmpBcType( IntField & bcTypeList );
};

class BasicRegion
{
public:
    int zid;
    int dir;         //�߽��淽��:0, 1, 2��Ӧ��i, j, k
    int outerNormal; //���ұ߽�-1, 1��Ӧ�����ұ߽�
    int start[ 3 ], end[ 3 ];
    int lr[ 3 ];     //���ұ߽�-1, 1��Ӧ�����ұ߽�
public:
    void SetRegion( int ist, int ied, int jst, int jed, int kst, int ked );
};

class BcRegion
{
public:
    BcRegion( int zid, int rid );
    ~BcRegion();
public:
    int rid;                         //region id
    int    bcType;                         //boundary type
    string regionName;               //boundary name
public:
    BasicRegion * s;
    BasicRegion * t;
public:
    void GetNormalizeIJKRegion( int & ist, int & ied, int & jst, int & jed, int & kst, int & ked );
    int ComputeRegionCells();
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
};

EndNameSpace