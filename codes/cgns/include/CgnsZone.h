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
#include "HXArray.h"
#include "GridDef.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsBase;
class CgnsCoor;
class NodeMesh;
class CgnsZsection;
class CgnsZbc;
class GridElem;
class PointFactory;
class ElemFeature;
class FaceSolver;
class CgnsBase;

class CgnsZone
{
public:
    CgnsZone( CgnsBase * cgnsBase );
    ~CgnsZone();
public:
    CgnsBase * cgnsBase;
    CgnsCoor * cgnsCoor;
    CgnsZsection * cgnsZsection;
    CgnsZbc * cgnsZbc;
public:
    string zoneName;
    int zId;

    ZoneType_t cgnsZoneType;

    int volBcType;

    IntField l2g;

    CgInt isize[ 9 ];
public:
    void CopyISize( CgInt * isize );
    void SetVolBcType( int volBcType );
    int GetVolBcType();
public:
    void Create();
    void SetPeriodicBc();
    void InitElement( GridElem * ge );
    void ConstructCgnsGridPoints( PointFactory * point_factory );
    void SetElementTypeAndNode( ElemFeature  * elem_feature );
    void InitLgMapping();
    void ConvertToInnerDataStandard();
public:
    void ScanBcFace( FaceSolver * face_solver );
    void GetElementNodeId( CgInt eId, CgIntField & eNodeId );
    void ReadCgnsGrid();
    void ReadCgnsZoneAttribute();
    void ReadCgnsZoneType();
    void ReadCgnsZoneNameAndGeneralizedDimension();
    void SetDimension();
    void ReadElementConnectivities();
    void ReadNumberOfCgnsSections();
    void CreateCgnsSections();
    void ReadCgnsSections();
    void ReadCgnsGridCoordinates();
    void ReadCgnsGridBoundary();
    void ProcessPeriodicBc();
    void ReadCgnsZoneBasicInfo();
public:
    void SetElemPosition();
public:
    CgInt GetNI() const;
    CgInt GetNJ() const;
    CgInt GetNK() const;
public:
    bool ExistSection( const string & sectionName );
    void GoToZone();
    void GoToNode( const string & nodeName, int ith );
    void GoToNode( const string & nodeNamei, int ith, const string & nodeNamej, int jth );
    void GoToNode( const string & nodeNamei, int ith, const string & nodeNamej, int jth, const string & nodeNamek, int kth );
public:
    void ReadFlowEqn();
};

#endif

EndNameSpace