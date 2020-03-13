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

#include "CgnsBase.h"
#include "CgnsBaseUtil.h"
#include "CgnsZone.h"
#include "CgnsZoneUtil.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "CgnsFamilyBc.h"
#include "GridMediator.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

void ReadAllCgnsZones( CgnsBase * myCgnsBase, CgnsBase * cgnsBaseIn )
{
    cout << "** Reading CGNS Grid In Base " << myCgnsBase->baseId << "\n";
    cout << "   numberOfCgnsZones       = " << myCgnsBase->nZones << "\n\n";

    for ( int iZone = 0; iZone < myCgnsBase->nZones; ++ iZone )
    {
        cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << myCgnsBase->nZones << "\n";
        CgnsZone * cgnsZone = myCgnsBase->GetCgnsZone( iZone );
        CgnsZone * cgnsZoneIn = cgnsBaseIn->GetCgnsZone( iZone );
        ONEFLOW::ReadCgnsGrid( cgnsZone, cgnsZoneIn );
    }
}

void ReadNumberOfCgnsZones( CgnsBase * myCgnsBase, CgnsBase * cgnsBaseIn )
{
    myCgnsBase->nZones = cgnsBaseIn->nZones;
}

void ReadCgnsBaseBasicInfo( CgnsBase * myCgnsBase, CgnsBase * cgnsBaseIn )
{
    myCgnsBase->baseName = cgnsBaseIn->baseName;
    myCgnsBase->celldim  = cgnsBaseIn->celldim;
    myCgnsBase->phydim   = cgnsBaseIn->phydim;
}

void DumpBase( CgnsBase * myCgnsBase, GridMediator * gridMediator )
{
    GlobalGrid::SetCurrentGridMediator( gridMediator );

    myCgnsBase->DumpCgnsBaseBasicInfo();

    cout << " nZones = " << myCgnsBase->nZones << "\n";

    for ( int iZone = 0; iZone < myCgnsBase->nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = myCgnsBase->GetCgnsZone( iZone );
        Grid * grid = gridMediator->gridVector[ iZone ];
        ONEFLOW::DumpCgnsZone( cgnsZone, grid );
    }
}

void PrepareCgnsZone( CgnsBase * myCgnsBase, GridMediator * gridMediator )
{
    GlobalGrid::SetCurrentGridMediator( gridMediator );

    cout << " nZones = " << myCgnsBase->nZones << "\n";

    for ( int iZone = 0; iZone < myCgnsBase->nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = myCgnsBase->GetCgnsZone( iZone );
        Grid * grid = gridMediator->gridVector[ iZone ];
        ONEFLOW::PrepareCgnsZone( cgnsZone, grid );
    }
}

#endif
EndNameSpace