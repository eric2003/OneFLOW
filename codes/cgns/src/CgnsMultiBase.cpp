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
    volBcType = -1;
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

void CgnsMultiBase::ReadCgnsBaseBasicInfo()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->ReadCgnsBaseBasicInfo();
    }
}

void CgnsMultiBase::ReadNumberOfCgnsZones()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->ReadNumberOfCgnsZones();
    }
}

void CgnsMultiBase::ReadCgnsBaseBasicInfo( CgnsMultiBase * strCgnsMultiBase )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->GetCgnsBase( bId );

        cgnsBase->ReadCgnsBaseBasicInfo( cgnsBaseIn );
    }
}

void CgnsMultiBase::ReadNumberOfCgnsZones( CgnsMultiBase * strCgnsMultiBase )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->GetCgnsBase( bId );

        cgnsBase->ReadNumberOfCgnsZones( cgnsBaseIn );
    }
}

void CgnsMultiBase::ReadCgnsMultiBase()
{
    this->ReadNumCgnsBase();

    this->AllocateCgnsBase();

    this->InitCgnsBase();

    this->ReadCgnsBaseBasicInfo();

    this->ReadNumberOfCgnsZones();

    this->ComputeNumberOfTotalZones();

    this->AllocateCgnsZonesInEachCgnsBase();

    this->InitAllCgnsZonesInEachCgnsBase();

    this->ReadAllCgnsZonesInEachCgnsBase();
}

void CgnsMultiBase::DumpCgnsMultiBase( GridMediator * gridMediator )
{
    int nZone = gridMediator->numberOfZones;
    this->Create( nZone );

    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->DumpBase( gridMediator );
    }
}

void CgnsMultiBase::ReadCgnsMultiBase( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadNumCgnsBase( strCgnsMultiBase );

    this->AllocateCgnsBase();

    this->InitCgnsBase();

    this->ReadCgnsBaseBasicInfo( strCgnsMultiBase );
    this->ReadNumberOfCgnsZones( strCgnsMultiBase );

    this->ComputeNumberOfTotalZones();

    this->AllocateCgnsZonesInEachCgnsBase();

    this->InitAllCgnsZonesInEachCgnsBase();

    this->ReadAllCgnsZonesInEachCgnsBase( strCgnsMultiBase );

}

void CgnsMultiBase::ReadAllCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->ReadFamilySpecifiedBc();

        cgnsBase->ReadAllCgnsZones();
    }
}

void CgnsMultiBase::ReadAllCgnsZonesInEachCgnsBase( CgnsMultiBase * strCgnsMultiBase )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->GetCgnsBase( bId );

        cgnsBase->ReadAllCgnsZones( cgnsBaseIn );
    }
}

void CgnsMultiBase::Create( int nZones )
{
    this->InitDefaultCgnsBase();

    this->AllocateCgnsBase();

    this->InitCgnsBase();

    this->SetGeneralizedCgnsZoneScale( nZones );

    this->ComputeNumberOfTotalZones();

    this->AllocateCgnsZonesInEachCgnsBase();

    this->InitAllCgnsZonesInEachCgnsBase();
}

void CgnsMultiBase::InitDefaultCgnsBase()
{
    this->fileId = 1;
    this->nBases = 1;
}

void CgnsMultiBase::AllocateCgnsBase()
{
    zid1.resize( this->nBases );
    zid2.resize( this->nBases );

    baseVector.resize( this->nBases );

    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = new CgnsBase();

        baseVector[ bId - 1 ] = cgnsBase;
    }
}

void CgnsMultiBase::InitCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        cgnsBase->fileId = this->fileId;
        cgnsBase->baseId = bId;
    }
}

void CgnsMultiBase::SetGeneralizedCgnsZoneScale( int nZones )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->SetDefaultCgnsBaseBasicInformation();

        cgnsBase->nZones = nZones;
    }
}

void CgnsMultiBase::ComputeNumberOfTotalZones()
{
    this->nTZones = 0;

    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        cgnsBase->zst = this->nTZones;
        cgnsBase->zed = this->nTZones + cgnsBase->nZones - 1;
        int id = bId - 1;
        this->zid1[ id ] = cgnsBase->zst;
        this->zid2[ id ] = cgnsBase->zed;

        this->nTZones += cgnsBase->nZones;
    }
}

void CgnsMultiBase::AllocateCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->AllocateAllCgnsZonesInCurrentCgnsBase();
    }
}

void CgnsMultiBase::InitAllCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );

        cgnsBase->InitAllCgnsZonesInCurrentCgnsBase();
    }
}

void CgnsMultiBase::ConvertStrCgns2UnsCgnsGrid( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadCgnsMultiBase( strCgnsMultiBase );
}

CgnsBase * CgnsMultiBase::GetCgnsBase( int baseId )
{
    return baseVector[ baseId - 1 ];
}

CgnsZone * CgnsMultiBase::GetCgnsZone( int globalZoneId )
{
    int bId = this->FindBaseId( globalZoneId );

    //int zst = zid1[ bId - 1 ];

    CgnsBase * cgnsBase = this->GetCgnsBase( bId );
    int zst = cgnsBase->zst;
    int localZoneId = globalZoneId - zst;
    CgnsZone * cgnsZone = cgnsBase->cgnsZones[ localZoneId ]; 

    return cgnsZone;
}

int CgnsMultiBase::FindBaseId( int iZone )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( bId );
        //int zst = zid1[ bId - 1 ];
        //int z2 = zid2[ bId - 1 ];

        int zst = cgnsBase->zst;
        int zed = cgnsBase->zed;

        if ( ( zst <= iZone ) && ( iZone <= zed ) )
        {
            return bId;
        }
    }

    return - 1;
}
#endif
EndNameSpace