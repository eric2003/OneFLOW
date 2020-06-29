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
#include <fstream>
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class ElemFeature;

class CgnsSection
{
public:
    CgnsSection( CgnsZone * cgnsZone );
    ~CgnsSection();
public:
    int eType;
    CgInt startId, endId;
    int connSize;
    CgIntField connList;

    IntField eTypeList;
    IntField ePosList;
    int pos_shift;

    string sectionName;
    CgnsZone * cgnsZone;
    int id;
    int nCoor;
    int nElement;
    CgInt elementDataSize;

    int nbndry;
    int iparentflag;
    CgIntField iparentdata;
public:
    void ConvertToInnerDataStandard();
    void SetElementTypeAndNode( ElemFeature * elem_feature );
    CgInt * GetAddress( CgInt eId );
    void GetElementNodeId( CgInt eId, CgIntField & eNodeId );
public:
    void ReadCgnsSection();
    void ReadCgnsSectionInfo();
    void SetSectionInfo( const string & sectionName, int elemType, int startId, int endId );
    void CreateConnList();
    void CalcNumberOfSectionElements();
    void CalcCapacityOfCgnsConnectionList();
    void AllocateCgnsConnectionList();
    void ReadCgnsSectionConnectionList();
    void SetElemPosition();
    void SetElemPositionOri();
    void SetElemPositionMixed();
    bool IsMixedSection();
};

#endif

EndNameSpace