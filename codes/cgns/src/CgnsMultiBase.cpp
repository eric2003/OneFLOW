/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "CgnsMultiBase.h"
#include "CgnsBase.h"
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


CgnsMultiBase::CgnsMultiBase()
{
}

CgnsMultiBase::~CgnsMultiBase()
{
    ;
}

int CgnsMultiBase::GetSystemZoneType()
{
    IntSet zoneTypeSet;

    for ( int iZone = 0; iZone < nTZones; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        int zoneType = cgnsZone->cgnsZoneType;
        zoneTypeSet.insert( zoneType );
    }

    if ( zoneTypeSet.size() == 1 )
    {
        return * zoneTypeSet.begin();
    }
    return ZoneTypeUserDefined;
}


void CgnsMultiBase::ReadCgnsGrid()
{
    this->ReadCgnsGrid( grid_para.gridFile );
}

void CgnsMultiBase::DumpCgnsGrid( GridMediator * gridMediator )
{
    string fileName = gridMediator->targetFile;
    this->OpenCgnsFile( fileName, CG_MODE_WRITE );
    this->DumpCgnsMultiBase( gridMediator );
    this->CloseCgnsFile();
}

void CgnsMultiBase::ReadCgnsGrid( const string & fileName )
{
    this->OpenCgnsFile( fileName, CG_MODE_READ );

    this->ReadCgnsMultiBase();

    this->CloseCgnsFile();
}

void CgnsMultiBase::OpenCgnsFile( const string & fileName, int cgnsOpenMode )
{
    //Open the CGNS for reading and check if the file was found.
    string prjFileName = GetPrjFileName( fileName );
    if ( cg_open( prjFileName.c_str(), cgnsOpenMode, & this->fileId ) != CG_OK )
    {
        Stop( cg_get_error() );
    }
}

void CgnsMultiBase::CloseCgnsFile()
{
    //close CGNS file
    cg_close( this->fileId );
}


void CgnsMultiBase::ReadNumCgnsBase()
{
    //Determine the of bases in the grid
    cg_nbases( this->fileId, & this->nBases );
    cout << "   Total number of CGNS Base = " << this->nBases << "\n";
}

void CgnsMultiBase::ReadNumCgnsBase( CgnsMultiBase * strCgnsMultiBase )
{
    this->fileId = strCgnsMultiBase->fileId;
    this->nBases = strCgnsMultiBase->nBases;
}

void CgnsMultiBase::ReadCgnsMultiBase()
{
    this->ReadNumCgnsBase();

    this->InitCgnsBase();

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );

        cgnsBase->ReadCgnsBaseBasicInfo();
        cgnsBase->ReadNumberOfCgnsZones();
        cgnsBase->AllocateAllCgnsZonesInCurrentCgnsBase();
        cgnsBase->InitAllCgnsZonesInCurrentCgnsBase();
        cgnsBase->ReadFamilySpecifiedBc();
        cgnsBase->ReadAllCgnsZones();
    }

    this->ComputeNumberOfTotalZones();
}


void CgnsMultiBase::DumpCgnsMultiBase( GridMediator * gridMediator )
{
    int nZone = gridMediator->numberOfZones;
    this->CreateDefaultCgnsZones( nZone );

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );

        cgnsBase->DumpBase( gridMediator );
    }
}

void CgnsMultiBase::ReadCgnsMultiBase( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadNumCgnsBase( strCgnsMultiBase );


    this->InitCgnsBase();

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->GetCgnsBase( iBase );

        cgnsBase->ReadCgnsBaseBasicInfo( cgnsBaseIn );
        cgnsBase->ReadNumberOfCgnsZones( cgnsBaseIn );
        cgnsBase->AllocateAllCgnsZonesInCurrentCgnsBase();
        cgnsBase->InitAllCgnsZonesInCurrentCgnsBase();
        cgnsBase->ReadAllCgnsZones( cgnsBaseIn );

    }

    this->ComputeNumberOfTotalZones();

}

void CgnsMultiBase::CreateDefaultCgnsZones( int nZones )
{
    this->InitDefaultCgnsBase();

    this->InitCgnsBase();

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );

        cgnsBase->SetDefaultCgnsBaseBasicInfo();
        cgnsBase->nZones = nZones;

        cgnsBase->AllocateAllCgnsZonesInCurrentCgnsBase();
        cgnsBase->InitAllCgnsZonesInCurrentCgnsBase();
    }

    this->ComputeNumberOfTotalZones();
}

void CgnsMultiBase::InitDefaultCgnsBase()
{
    this->fileId = 1;
    this->nBases = 1;
}

void CgnsMultiBase::InitCgnsBase()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = new CgnsBase();
        baseVector.push_back( cgnsBase );

        cgnsBase->fileId = this->fileId;
        cgnsBase->baseId = iBase + 1;
    }
}

void CgnsMultiBase::ComputeNumberOfTotalZones()
{
    this->nTZones = 0;

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        cgnsBase->zst = this->nTZones;
        cgnsBase->zed = this->nTZones + cgnsBase->nZones - 1;

        this->nTZones += cgnsBase->nZones;
    }
}

void CgnsMultiBase::ConvertStrCgns2UnsCgnsGrid( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadCgnsMultiBase( strCgnsMultiBase );
}

CgnsBase * CgnsMultiBase::GetCgnsBase( int iBase )
{
    return baseVector[ iBase ];
}

CgnsZone * CgnsMultiBase::GetCgnsZone( int globalZoneId )
{
    CgnsZone * cgnsZone = this->FindGlobalCgnsZone( globalZoneId );
    return cgnsZone;
}

CgnsZone * CgnsMultiBase::FindGlobalCgnsZone( int globalZoneId )
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        CgnsZone * cgnsZone = cgnsBase->FindGlobalCgnsZone( globalZoneId );
        if ( cgnsZone )
        {
            return cgnsZone;
        }
    }

    return 0;
}

#endif
EndNameSpace