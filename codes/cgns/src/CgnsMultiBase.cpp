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
        CgnsZone * cgnsZone = this->GetZone( iZone );
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

void CgnsMultiBase::ReadGeneralizedCgnsZoneScale()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = baseVector[ bId - 1 ];

        cgnsBase->ReadCgnsBaseBasicInfo();

        cgnsBase->ReadNumberOfCgnsZones();
    }
}

void CgnsMultiBase::ReadGeneralizedCgnsZoneScale( CgnsMultiBase * strCgnsMultiBase )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = this->baseVector[ bId - 1 ];
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->baseVector[ bId - 1 ];

        cgnsBase->ReadCgnsBaseBasicInfo( cgnsBaseIn );

        cgnsBase->ReadNumberOfCgnsZones( cgnsBaseIn );
    }
}

void CgnsMultiBase::ReadCgnsMultiBase()
{
    this->ReadNumCgnsBase();

    this->AllocateCgnsBase();

    this->InitCgnsBase();

    this->ReadGeneralizedCgnsZoneScale();

    this->ComputeNumberOfTotalZones();

    this->AllocateCgnsZonesInEachCgnsBase();

    this->InitAllCgnsZonesInEachCgnsBase();

    this->ReadAllCgnsZonesInEachCgnsBase();
}

void CgnsMultiBase::ReadCgnsMultiBase( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadNumCgnsBase( strCgnsMultiBase );

    this->AllocateCgnsBase();

    this->InitCgnsBase();

    this->ReadGeneralizedCgnsZoneScale( strCgnsMultiBase );

    this->ComputeNumberOfTotalZones();

    this->AllocateCgnsZonesInEachCgnsBase();

    this->InitAllCgnsZonesInEachCgnsBase();

    this->ReadAllCgnsZonesInEachCgnsBase( strCgnsMultiBase );

}

void CgnsMultiBase::ReadAllCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int id = bId - 1;

        CgnsBase * cgnsBase = baseVector[ id ];

        cgnsBase->ReadFamilySpecifiedBc();

        cgnsBase->ReadAllCgnsZones();
    }
}

void CgnsMultiBase::ReadAllCgnsZonesInEachCgnsBase( CgnsMultiBase * strCgnsMultiBase )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int id = bId - 1;

        CgnsBase * cgnsBase = baseVector[ id ];
        CgnsBase * cgnsBaseIn = strCgnsMultiBase->baseVector[ id ];

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
        CgnsBase * cgnsBase = baseVector[ bId - 1 ];
        cgnsBase->fileId = this->fileId;
        cgnsBase->baseId = bId;
    }
}

void CgnsMultiBase::SetGeneralizedCgnsZoneScale( int nZones )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        CgnsBase * cgnsBase = baseVector[ bId - 1 ];

        cgnsBase->SetDefaultCgnsBaseBasicInformation();

        cgnsBase->nZones = nZones;
    }
}

void CgnsMultiBase::ComputeNumberOfTotalZones()
{
    this->nTZones = 0;

    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int id = bId - 1;
        CgnsBase * cgnsBase = baseVector[ id ];

        this->zid1[ id ] = this->nTZones;
        this->zid2[ id ] = this->nTZones + cgnsBase->nZones - 1;

        this->nTZones += cgnsBase->nZones;
    }
}

void CgnsMultiBase::AllocateCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int id = bId - 1;

        CgnsBase * cgnsBase = baseVector[ id ];

        cgnsBase->AllocateAllCgnsZonesInCurrentCgnsBase();
    }
}

void CgnsMultiBase::InitAllCgnsZonesInEachCgnsBase()
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int id = bId - 1;

        CgnsBase * cgnsBase = baseVector[ id ];

        cgnsBase->InitAllCgnsZonesInCurrentCgnsBase();
    }
}

void CgnsMultiBase::ConvertStrCgns2UnsCgnsGrid( CgnsMultiBase * strCgnsMultiBase )
{
    this->ReadCgnsMultiBase( strCgnsMultiBase );
}

CgnsZone * CgnsMultiBase::GetZone( int iZone )
{
    int bId = this->FindBaseId( iZone );

    int z1 = zid1[ bId - 1 ];

    CgnsZone * cgnsZone = baseVector[ bId - 1 ]->cgnsZones[ iZone - z1 ]; 

    return cgnsZone;
}

int CgnsMultiBase::FindBaseId( int iZone )
{
    for ( int bId = 1; bId <= this->nBases; ++ bId )
    {
        int z1 = zid1[ bId - 1 ];
        int z2 = zid2[ bId - 1 ];

        if ( ( z1 <= iZone ) && ( iZone <= z2 ) )
        {
            return bId;
        }
    }

    return - 1;
}
#endif
EndNameSpace