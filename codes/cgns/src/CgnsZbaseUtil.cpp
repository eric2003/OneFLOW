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

#include "CgnsZbase.h"
#include "CgnsZbaseUtil.h"
#include "CgnsBase.h"
#include "CgnsBaseUtil.h"
#include "CgnsZone.h"
#include "StrUtil.h"
#include "Stop.h"
#include "Prj.h"
#include "Dimension.h"
#include "GridPara.h"
#include "GridMediator.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS


void ReadNumCgnsBase( CgnsZbase * myCgnsZbase, CgnsZbase * strCgnsMultiBase )
{
    myCgnsZbase->fileId = strCgnsMultiBase->fileId;
    myCgnsZbase->nBases = strCgnsMultiBase->nBases;
}

void ReadCgnsMultiBase( CgnsZbase * myCgnsZbase, CgnsZbase * strCgnsMultiBase )
{
    ONEFLOW::ReadNumCgnsBase( myCgnsZbase, strCgnsMultiBase );

    myCgnsZbase->InitCgnsBase();

    for ( int iBase = 0; iBase < myCgnsZbase->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = myCgnsZbase->GetCgnsBase( iBase );
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->GetCgnsBase( iBase );

        ONEFLOW::ReadCgnsBaseBasicInfo( cgnsBase, cgnsBaseIn );
        ONEFLOW::ReadNumberOfCgnsZones( cgnsBase, cgnsBaseIn );
        cgnsBase->AllocateAllCgnsZones();
        ONEFLOW::ReadAllCgnsZones( cgnsBase, cgnsBaseIn );

    }
}

void ConvertStrCgns2UnsCgnsGrid( CgnsZbase * myCgnsZbase, CgnsZbase * strCgnsMultiBase )
{
    ONEFLOW::ReadCgnsMultiBase( myCgnsZbase, strCgnsMultiBase );
}

void CreateDefaultCgnsZones( CgnsZbase * myCgnsZbase, GridMediatorS * gridMediatorS )
{
    myCgnsZbase->fileId = 1;
    myCgnsZbase->nBases = gridMediatorS->GetSize();

    myCgnsZbase->InitCgnsBase();

    for ( int iBase = 0; iBase < myCgnsZbase->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = myCgnsZbase->GetCgnsBase( iBase );
        GridMediator * gridMediator = gridMediatorS->GetGridMediator( iBase );

        cgnsBase->SetDefaultCgnsBaseBasicInfo();
        cgnsBase->nZones = gridMediator->numberOfZones;

        cgnsBase->AllocateAllCgnsZones();
    }
}

void DumpCgnsMultiBase( CgnsZbase * myCgnsZbase, GridMediatorS * gridMediatorS )
{
    ONEFLOW::CreateDefaultCgnsZones( myCgnsZbase, gridMediatorS );

    for ( int iBase = 0; iBase < myCgnsZbase->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = myCgnsZbase->GetCgnsBase( iBase );
        GridMediator * gridMediator = gridMediatorS->GetGridMediator( iBase );
        ONEFLOW::DumpBase( cgnsBase, gridMediator );
    }
}

void DumpCgnsGrid( CgnsZbase * myCgnsZbase, GridMediatorS * gridMediators )
{
    string fileName = gridMediators->GetTargetFile();
    myCgnsZbase->OpenCgnsFile( fileName, CG_MODE_WRITE );
    ONEFLOW::DumpCgnsMultiBase( myCgnsZbase, gridMediators );
    myCgnsZbase->CloseCgnsFile();
}

void PrepareCgnsZone( CgnsZbase * myCgnsZbase, GridMediatorS * gridMediatorS )
{
    for ( int iBase = 0; iBase < myCgnsZbase->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = myCgnsZbase->GetCgnsBase( iBase );
        GridMediator * gridMediator = gridMediatorS->GetGridMediator( iBase );

        ONEFLOW::PrepareCgnsZone( cgnsBase, gridMediator );
    }
}

#endif
EndNameSpace