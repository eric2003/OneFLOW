/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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


BeginNameSpace( ONEFLOW )

class GridMediator;
struct GridConfig;
class TextFileParser;
const int MAX_VTK_TYPE = 100;
class VTK_TYPE
{
public:
    VTK_TYPE();
    ~VTK_TYPE();
public:
    static int VERTEX       ;
    static int LINE         ;
    static int TRIANGLE     ;
    static int QUADRILATERAL;
    static int TETRAHEDRON  ;
    static int HEXAHEDRON   ;
    static int PRISM        ;
    static int PYRAMID      ;
};

class VTK_CgnsMap
{
public:
    VTK_CgnsMap();
    ~VTK_CgnsMap();
public:
    std::map< int, int > vtk2Cgns;
public:
    void Init();
};

extern VTK_CgnsMap vtk_CgnsMap;

class Marker
{
public:
    Marker() = default;
    ~Marker() = default;
public:
    std::string name;
    std::string bcName;
    int cgns_bcType{ 0 };
    LinkField elems;
    IntField eTypes;
    int nElem{ 0 };
};

class SecMarker
{
public:
    SecMarker() = default;
    ~SecMarker() = default;
public:
    int vtk_type{ 0 };
    int cgns_type{ 0 };
    std::string name;
    int nElem{ 0 };
    LinkField elems;
};

class SecMarkerManager
{
public:
    SecMarkerManager();
    ~SecMarkerManager() = default;
public:
    int nType;
    HXVector< SecMarker > data;
public:
    void Alloc( int nType );
    int CalcTotalElem();
};

class MarkerManager
{
public:
    MarkerManager();
    ~MarkerManager() = default;
public:
    int nMarker;
    HXVector< Marker > markerList;

    IntField types;
    LinkField l2g;
public:
    void CreateMarkerList( int nMarker );
    void CalcSecMarker( SecMarkerManager * secMarkerManager );
};

class Su2Grid;
class VolumeSecManager
{
public:
    VolumeSecManager();
    ~VolumeSecManager();
public:
    int nSection;

    IntField types;
    LinkField l2g;
public:
    void CalcVolSec( Su2Grid* su2Grid, SecMarkerManager * secMarkerManager );
};

class Su2Bc
{
public:
    Su2Bc();
    ~Su2Bc();
public:
    std::set< std::string > bcList;
    std::map<std::string, std::string> bcMap;
    std::map<std::string, int> bcNameToValueMap;
public:
    void Init();
    void AddBc( std::string &geoName, std::string &bcName);
    void Process(StringField& markerBCNameList, StringField& markerNameList);
    std::string GetBcName( std::string& geoName );
    int GetCgnsBcType(std::string& geoName);
};

class CgnsZone;

class Su2Grid
{
public:
    Su2Grid();
    ~Su2Grid();
public:
    void ReadSu2Grid( GridMediator * gridMediator );
    void ReadSu2GridAscii( const std::string & fileName, const std::string & caseDir );
    void Su2ToOneFlowGrid( const GridConfig & config, const std::string & caseDir );
    void MarkBoundary( std::string & su2cfgFile, const std::string & caseDir );
    void FillSU2CgnsZone( CgnsZone * cgnsZone );
public:
    int ndim;
    int nPoin, nElem;
    IntField vtkmap;
    IntField elemVTKType;
    IntField elemId;
    LinkField elems;
    RealField xN, yN, zN;
    MarkerManager mmark;
    VolumeSecManager volSec;
    Su2Bc su2Bc;
public:
    int nZone;

};

void Su2ToOneFlowGrid( Su2Grid & su2Grid );

EndNameSpace
