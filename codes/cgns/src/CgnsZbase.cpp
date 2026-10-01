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

#include "CgnsZbase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "CgnsFile.h"
#include "Fatal.h"
#include <iostream>
#include <utility>

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbase::CgnsZbase()
{
    this->nBases = 0;
    this->cgnsFile = std::make_unique< CgnsFile >();
}

CgnsZbase::~CgnsZbase() = default;

void CgnsZbase::FreeCgnsBases()
{
    this->baseVector.clear();
    this->nBases = 0;
}

void CgnsZbase::OpenCgnsFile( const std::string & fileName, int cgnsOpenMode )
{
    this->cgnsFile->OpenCgnsFile( fileName, cgnsOpenMode );
}

void CgnsZbase::CloseCgnsFile()
{
    this->cgnsFile->CloseCgnsFile();
}

int CgnsZbase::GetNZones()
{
    int nZones = 0;
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        nZones += cgnsBase->GetNZones();
    }
    return nZones;
}

void CgnsZbase::CreateCgnsZones( int nZones )
{
    CgnsBase * cgnsBase = this->CreateCgnsBase();
    cgnsBase->CreateCgnsZones( nZones );
}

CgnsZone * CgnsZbase::CreateCgnsZone()
{
    CgnsBase * cgnsBase = 0;
    if ( this->nBases == 0 )
    {
        cgnsBase = this->CreateCgnsBase();
    }
    else
    {
        cgnsBase = this->GetCgnsBase( 0 );
    }
    
    CgnsZone * cgnsZone = cgnsBase->CreateCgnsZone();
    return cgnsZone; 
}

int CgnsZbase::GetSystemZoneType()
{
    IntSet zoneTypeSet;
    int nTZones = this->GetNZones();
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

void CgnsZbase::ReadCgnsGrid( const std::string & fileName )
{
    this->OpenCgnsFile( fileName, CG_MODE_READ );
    this->ReadCgnsMultiBase();
    this->CloseCgnsFile();
}

void CgnsZbase::DumpCgnsMultiBase()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );

        cgnsBase->DumpCgnsBaseBasicInfo();
        cgnsBase->DumpAllCgnsZones();
    }
}


void CgnsZbase::ReadNumCgnsBase()
{
    //Determine the of bases in the grid
    cg_nbases( this->cgnsFile->fileId, & this->nBases );
    std::cout << "   Total number of CGNS Base = " << this->nBases << "\n";
}

void CgnsZbase::ConvertToInnerDataStandard()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        cgnsBase->ConvertToInnerDataStandard();
    }
}

void CgnsZbase::ProcessCgnsBases()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->GetCgnsBase( iBase );
        cgnsBase->ProcessCgnsZones();
    }
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
    this->AddCgnsBase( std::unique_ptr< CgnsBase >( cgnsBase ) );
}

void CgnsZbase::AddCgnsBase( std::unique_ptr< CgnsBase > cgnsBase )
{
    CgnsBase * base = cgnsBase.get();
    baseVector.push_back( std::move( cgnsBase ) );
    int baseId = baseVector.size();
    base->cgnsFile = this->cgnsFile.get();
    base->baseId = baseId;
}

CgnsBase * CgnsZbase::CreateCgnsBase()
{
    auto cgnsBase = std::make_unique< CgnsBase >( this->cgnsFile.get() );
    CgnsBase * base = cgnsBase.get();
    this->AddCgnsBase( std::move( cgnsBase ) );
    return base;
}

void CgnsZbase::InitCgnsBase()
{
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        this->CreateCgnsBase();
    }
}

CgnsBase * CgnsZbase::GetCgnsBase( int iBase )
{
    return baseVector[ iBase ].get();
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
