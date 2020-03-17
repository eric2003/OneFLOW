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
#include "CgnsZone.h"
#include "CgnsZoneUtil.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "CgnsFamilyBc.h"
#include "CgnsVariable.h"

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

CgnsBase::CgnsBase()
{
    this->familyBc = 0;
    this->freeFlag = false;
}

CgnsBase::~CgnsBase()
{
    delete this->familyBc;
    if ( this->freeFlag )
    {
        this->FreeZoneList();
    }
}

void CgnsBase::FreeZoneList()
{
    for ( int i = 0; i < cgnsZones.size(); ++ i )
    {
        delete cgnsZones[ i ];
    }
}


CgnsZone * CgnsBase::GetCgnsZone( int iZone )
{
    //iZone base on 0
    return this->cgnsZones[ iZone ];
}

CgnsZone * CgnsBase::GetCgnsZoneByName( const string & zoneName )
{
    map< string, int >::iterator iter;
    iter = zoneNameMap.find( zoneName );
    int iZone = iter->second - 1;
    return this->GetCgnsZone( iZone );
}

int CgnsBase::GetNZone()
{
    return this->cgnsZones.size();
}

void CgnsBase::SetDefaultCgnsBaseBasicInfo()
{
    //this->celldim = Dim::dimension;
    //this->phydim  = Dim::dimension;

    this->celldim = THREE_D;
    this->phydim  = THREE_D;
  
    this->baseName = ONEFLOW::AddString( "Base", this->baseId );
}

void CgnsBase::AddCgnsZone( CgnsZone * cgnsZone )
{
    cgnsZones.push_back( cgnsZone );
    int zId = cgnsZones.size();
    cgnsZone->zId = zId;
}

void CgnsBase::AllocateAllCgnsZones()
{
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = new CgnsZone( this );

        this->AddCgnsZone( cgnsZone );

        cgnsZone->Create();
    }
}

void CgnsBase::ReadCgnsBaseBasicInfo()
{
    CgnsTraits::char33 cgnsBaseName;

    //Check the cell and physical dimensions of the bases.
    cg_base_read( this->fileId, this->baseId, cgnsBaseName, & this->celldim, & this->phydim );
    this->baseName = cgnsBaseName;
    cout << "   baseId = " << this->baseId << " baseName = " << cgnsBaseName << "\n";
    cout << "   cell dim = " << this->celldim << " physical dim = " << this->phydim << "\n";
}

void CgnsBase::DumpCgnsBaseBasicInfo()
{
    cg_base_write( this->fileId, this->baseName.c_str(), this->celldim, this->phydim, &this->baseId );
    cout << " baseId = " << this->baseId << " baseName = " << this->baseName << "\n";
}

void CgnsBase::ReadNumberOfCgnsZones()
{
    //Read the number of zones in the grid.
    cg_nzones( this->fileId, this->baseId, & this->nZones );
}

void CgnsBase::ConstructZoneNameMap()
{
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        zoneNameMap[ cgnsZone->zoneName ] = cgnsZone->zId;
    }
}

void CgnsBase::ReadAllCgnsZones()
{
    cout << "** Reading CGNS Grid In Base " << this->baseId << "\n";
    cout << "   Reading CGNS Family Specified BC \n";
    this->ReadFamilySpecifiedBc();
    cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->ReadCgnsGrid();
    }

    this->ConstructZoneNameMap();

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        cout << "cgnsZone->SetPeriodicBc\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->SetPeriodicBc();
    }
}

void CgnsBase::SetFamilyBc( BCType_t & bcType, const string & bcRegionName )
{
    this->familyBc->SetFamilyBc( bcType, bcRegionName );
}

void CgnsBase::ReadFamilySpecifiedBc()
{
    this->familyBc = new CgnsFamilyBc( this );
    this->familyBc->ReadFamilySpecifiedBc();
}

void CgnsBase::WriteZoneInfo( const string & zoneName, ZoneType_t zoneType, cgsize_t * isize )
{
    int cgzone;
    cg_zone_write( this->fileId, this->baseId, zoneName.c_str(), isize, zoneType, & cgzone );
    this->freeFlag = true;

    CgnsZone * cgnsZone = new CgnsZone( this );
    cgnsZone->zoneName = zoneName;
    cgnsZone->cgnsZoneType = zoneType;
    cgnsZone->CopyISize( isize );
    this->cgnsZones.push_back( cgnsZone );    ;
}

void CgnsBase::GoToBase()
{
    cg_goto( this->fileId, this->baseId, "end" );
}

void CgnsBase::GoToNode( const string & nodeName, int ith )
{
    cg_goto( this->fileId, this->baseId, nodeName.c_str(), ith, NULL );
}

void CgnsBase::ReadArray()
{
    this->GoToBase();
    int nUserData = -1;
    cg_nuser_data( & nUserData );

    cout << " nUserData = " << nUserData << "\n";

    for ( int iData = 0; iData < nUserData; ++ iData )
    {
        int iDataId = iData + 1;
        this->GoToBase();

        char name[ 33 ];
        cg_user_data_read( iDataId, name );

        cout << " iData = " << iData << " user data name = " << name << "\n";

        this->GoToNode( "UserDefinedData_t", iDataId );
        int narrays = -1;
        cg_narrays( & narrays );

        cout << " narrays = " << narrays << "\n";

        CgnsZVector cgnsZVector;
        cgnsZVector.ReadArray( narrays );
    }
}


#endif
EndNameSpace