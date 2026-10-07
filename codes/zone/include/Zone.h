/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#pragma once
#include "HXDefine.h"
#include "GridHandles.h"
#include <fstream>
#include <memory>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

class Grid;
class UnsGrid;
class ScalarGrid;
class InterFaceTopo;

class Zone
{
public:
    Zone();
    ~Zone();
public:
    static int flag_test_grid;
    // One Grids list per zone id (multigrid levels stored as successive unique_ptrs).
    static std::vector< Grids > globalGrids;
    static int nLocalZones;
    static std::unique_ptr< InterFaceTopo > interfaceTopo;
    static void AddGrid( int zid, std::unique_ptr< Grid > grid );
    static void ReleaseGrids();
    static void Reset();
    static InterFaceTopo & GetInterfaceTopo();
    static void InitInterfaceTopo();
    static void InitLayout( StringField & fileNameList );
    static void InitLayout( StringField & fileNameList, const std::string & caseDir );
    static void ReadGrid( StringField & fileNameList );
    static void ReadGrid( StringField & fileNameList, const std::string & caseDir );
    static void NormalizeLayout();
public:
    static Grid * GetGrid( int zid, int gl = 0 );
    static Grid & GetGridReference( int zid, int gl = 0 );
    static Grid * GetGrid();
    static Grid & GetGridReference();
    static Grid * GetCGrid( Grid * grid );
    static Grid * GetFGrid( Grid * grid );
    static UnsGrid * GetUnsGrid();
public:
    static void AddScalarGrid( int zid, std::unique_ptr< ScalarGrid > grid );
    static ScalarGrid * GetScalarGrid( int iZone );
    static ScalarGrid * GetScalarGrid();
    static ScalarGrid & GetScalarGridReference( int iZone );
    static ScalarGrid & GetScalarGridReference();
public:
    static int GetNumberOfZoneNeighbors( int zoneId );
    static int GetNeighborZoneId( int zoneId, int iNeighbor );
};

EndNameSpace
