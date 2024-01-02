/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "GridDef.h"
#include <fstream>
#include <string>


BeginNameSpace( ONEFLOW )

class Grid;
class UnsGrid;
class ScalarGrid;

class Zone
{
public:
    Zone();
    ~Zone();
public:
    static int flag_test_grid;
    static HXVector< Grids * > globalGrids;
    static int nLocalZones;
    static void AddGrid( int zid, Grid * grid );
    static void InitLayout( StringField & fileNameList );
    static void ReadGrid( StringField & fileNameList );
    static void NormalizeLayout();
public:
    static Grid * GetGrid( int zid, int gl = 0 );
    static Grid * GetGrid();
    static Grid * GetCGrid( Grid * grid );
    static Grid * GetFGrid( Grid * grid );
    static UnsGrid * GetUnsGrid();
public:
    static void AddScalarGrid( int zid, ScalarGrid * grid );
    static ScalarGrid * GetScalarGrid( int iZone );
    static ScalarGrid * GetScalarGrid();
public:
    static int GetNumberOfZoneNeighbors( int zoneId );
    static int GetNeighborZoneId( int zoneId, int iNeighbor );
};

EndNameSpace
