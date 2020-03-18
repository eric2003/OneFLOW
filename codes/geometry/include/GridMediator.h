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
#include "GridDef.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class GridMediator
{
public:
    GridMediator ();
    ~GridMediator();
public:
    Grids gridVector;
    int numberOfZones;
    int readGridType;
    string gridFile;
    string bcFile;
    string targetFile;
    string gridType;
public:
    void ReadGrid();
    void ReadGridgen();
    void ReadPlot3D();
    void ReadPlot3DCoor();
public:
    void AddDefaultName();
};

class ZgridMediator
{
public:
    ZgridMediator();
    ~ZgridMediator();
public:
    HXVector< GridMediator * > gm;
    bool flag;
public:
    void AddGridMediator( GridMediator * gridMediator );
    GridMediator * GetGridMediator( int iGridMediator );
    int GetSize();
    string GetTargetFile();
    void SetDeleteFlag( bool flag );
public:
    void CreateSimple( int nZone );
    void ReadGrid();
};

class GlobalGrid
{
public:
    GlobalGrid();
    ~GlobalGrid();
public:
    static GridMediator * gridMediator;
    static Grid * GetGrid( int zoneId );
    static void SetCurrentGridMediator( GridMediator * gridMediatorIn );
};

EndNameSpace