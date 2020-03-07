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

class Grid;
class StrGrid;
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
    void ReadCgnsGrid( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneAttribute();
    void DumpCgnsZoneAttribute( Grid * grid );
    void ReadCgnsZoneAttribute( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneType();
    void DumpCgnsZoneType( Grid * grid );
    void ReadCgnsZoneType( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneNameAndGeneralizedDimension();
    void DumpCgnsZoneNameAndGeneralizedDimension( Grid * gridIn );
    void ReadCgnsZoneNameAndGeneralizedDimension( CgnsZone * cgnsZoneIn );
    void SetDimension();
    void SetDimension( CgnsZone * cgnsZoneIn );
    void ReadElementConnectivities();
    void ReadElementConnectivities( CgnsZone * cgnsZoneIn );
    void ReadNumberOfCgnsSections();
    void ReadNumberOfCgnsSections( CgnsZone * cgnsZoneIn );
    void CreateCgnsSections();
    void ReadCgnsSections();
    void ReadCgnsGridCoordinates();
    void DumpCgnsGridCoordinates( Grid * grid );
    void ReadCgnsGridCoordinates( CgnsZone * cgnsZoneIn );
    void ReadCgnsGridBoundary();
    void DumpCgnsGridBoundary( Grid * grid );
    void ProcessPeriodicBc();
    void DumpCgnsZone( Grid * grid );
    void FillISize( Grid * gridIn );
    void FillISize( int ni, int nj, int nk, int dimension );
    void PrepareCgnsZone( Grid * grid );
public:
    void SetElemPosition();
public:
    CgInt GetNI() const;
    CgInt GetNJ() const;
    CgInt GetNK() const;
public:
    //void GetStrZonePara( int & s1, int & e1, int & s2, int & e2, int & etype1, int & etype2 );
public:
    bool ExistSection( const string & sectionName );
};

#endif

EndNameSpace