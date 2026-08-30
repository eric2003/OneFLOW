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

#include "GridUtils.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "UnsGrid.h"
#include "FaceTopo.h"


BeginNameSpace( ONEFLOW )

int GetNumberOfSolidCells( UnsGrid * grid )
{
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
    bcRecord->CreateBcTypeRegion();

    BcInfo * bcInfo = bcRecord->bcInfo;

    int nRegion = bcInfo->bcType.size();

    int nSolidCells = 0;
    for ( int ir = 0; ir < nRegion; ++ ir )
    {
        int bcType = bcInfo->bcType[ ir ];
        if ( bcType != BC::SOLID_SURFACE ) continue;

        int nBCFace = bcInfo->bcFace[ ir ].size();
		nSolidCells += nBCFace;
    }

    return nSolidCells;
}

EndNameSpace
