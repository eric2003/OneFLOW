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
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "Stop.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS


CgnsZbase::CgnsZbase()
{
    this->fileId = 1;
    this->nBases = 0;
}

CgnsZbase::~CgnsZbase()
{
    ;
}

int CgnsZbase::GetNZone()
{
    int nZone = 0;
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        nZone += cgnsBase->GetNZone();
    }
    return nZone;
}

int CgnsZbase::GetSystemZoneType()
{
    IntSet zoneTypeSet;
    int nTZones = this->GetNZone();
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

void CgnsZbase::ReadCgnsGrid( const string & fileName )
{
    this->OpenCgnsFile( fileName, CG_MODE_READ );

    this->ReadCgnsMultiBase();

    this->CloseCgnsFile();
}

void CgnsZbase::OpenCgnsFile( const string & fileName, int cgnsOpenMode )
{
    //Open the CGNS for reading and check if the file was found.
    if ( cg_open( fileName.c_str(), cgnsOpenMode, & this->fileId ) != CG_OK )
    {
        Stop( cg_get_error() );
    }
}

void CgnsZbase::CloseCgnsFile()
{
    //close CGNS file
    cg_close( this->fileId );
}


void CgnsZbase::ReadNumCgnsBase()
{
    //Determine the of bases in the grid
    cg_nbases( this->fileId, & this->nBases );
    cout << "   Total number of CGNS Base = " << this->nBases << "\n";
}

void CgnsZbase::ReadCgnsMultiBase()
{
    this->ReadNumCgnsBase();

    this->InitCgnsBase();

    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );

        cgnsBase->ReadCgnsBaseBasicInfo();
        cgnsBase->ReadNumberOfCgnsZones();
        cgnsBase->AllocateAllCgnsZones();
        cgnsBase->ReadAllCgnsZones();
    }
}

void CgnsZbase::AddCgnsBase( CgnsBase * cgnsBase )
{
    baseVector.push_back( cgnsBase );
    int baseId = baseVector.size();
    cgnsBase->fileId = this->fileId;
    cgnsBase->baseId = baseId;
}

void CgnsZbase::InitCgnsBase()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = new CgnsBase();
        this->AddCgnsBase( cgnsBase );
    }
}

CgnsBase * CgnsZbase::GetCgnsBase( int iBase )
{
    return baseVector[ iBase ];
}

CgnsZone * CgnsZbase::GetCgnsZone( int globalZoneId )
{
    CgnsZone * cgnsZone = this->GetMultiBaseCgnsZone( 0, globalZoneId );
    return cgnsZone;
}

CgnsZone * CgnsZbase::GetMultiBaseCgnsZone( int iBase, int iZone )
{
    CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
    CgnsZone * cgnsZone = cgnsBase->GetCgnsZone( iZone );
    return cgnsZone;
}

#endif
EndNameSpace