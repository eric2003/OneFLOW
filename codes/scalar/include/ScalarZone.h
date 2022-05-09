/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

class ScalarGrid;

class ScalarZone
{
public:
    ScalarZone();
    ~ScalarZone();
public:
    static int nLocalZones;
    static HXVector< ScalarGrid * > scalar_grids;
public:
    static void Allocate();
    static void DeAllocate();
    static void AddGrid( int zid, ScalarGrid * grid );
    static ScalarGrid * GetGrid( int iZone );
    static ScalarGrid * GetGrid();

};

EndNameSpace
