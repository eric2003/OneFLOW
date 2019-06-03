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

#include "CgnsBase.h"
#include "CgnsZone.h"
#include "BasicIO.h"
#include "Dimension.h"

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

CgnsBase::CgnsBase()
{
}

CgnsBase::~CgnsBase()
{
}

CgnsZone * CgnsBase::GetCgnsZone( int zoneId )
{
	int id = zoneId - 1;
	return this->cgnsZones[ id ];
}

CgnsZone * CgnsBase::GetCgnsZone( const string & zoneName )
{
	map< string, int >::iterator iter;
	iter = zoneNameMap.find( zoneName );
	int zoneId = iter->second;
	return this->GetCgnsZone( zoneId );
}

void CgnsBase::SetDefaultCgnsBaseBasicInformation()
{
    this->celldim = Dim::dimension;
	this->phydim  = Dim::dimension;

	this->baseName = ONEFLOW::AddString( "Base", this->baseId );
}

void CgnsBase::AllocateAllCgnsZonesInCurrentCgnsBase()
{
	cgnsZones.resize( nZones );

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
        CgnsZone * cgnsZone = new CgnsZone( this );

		cgnsZones[ iZone ] = cgnsZone;

        cgnsZone->Create();
	}
}

void CgnsBase::InitAllCgnsZonesInCurrentCgnsBase()
{
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
        CgnsZone * cgnsZone = cgnsZones[ iZone ];
		cgnsZone->zId = iZone + 1;
	}
}

void CgnsBase::ReadCgnsBaseBasicInfo()
{
	CgnsTraits::char33 cgnsBaseName;

	//Check the cell and physical dimensions of the bases.
	cg_base_read( this->fileId, this->baseId, cgnsBaseName, & this->celldim, & this->phydim );
    this->baseName = cgnsBaseName;
	cout << "   baseId = " << this->baseId << " baseName = " << cgnsBaseName << "\n";
}

void CgnsBase::ReadCgnsBaseBasicInfo( CgnsBase * cgnsBaseIn )
{
    this->baseName = cgnsBaseIn->baseName;
    this->celldim  = cgnsBaseIn->celldim;
    this->phydim   = cgnsBaseIn->phydim;
}

void CgnsBase::ReadNumberOfCgnsZones()
{
	//Read the number of zones in the grid.
	cg_nzones( this->fileId, this->baseId, & this->nZones );
}

void CgnsBase::ReadNumberOfCgnsZones( CgnsBase * cgnsBaseIn )
{
	this->nZones = cgnsBaseIn->nZones;
}

void CgnsBase::ConstructZoneNameMap()
{
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		CgnsZone * cgnsZone = this->cgnsZones[ iZone ];
		zoneNameMap[ cgnsZone->zoneName ] = cgnsZone->zId;
	}
}

void CgnsBase::ReadAllCgnsZones()
{
    cout << "** Reading CGNS Grid In Base " << this->baseId << "\n";
    cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
        cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
		CgnsZone * cgnsZone = this->cgnsZones[ iZone ];
		cgnsZone->ReadCgnsGrid();
	}

	this->ConstructZoneNameMap();

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
		cout << "cgnsZone->SetPeriodicBc\n";
		CgnsZone * cgnsZone = this->cgnsZones[ iZone ];
		cgnsZone->SetPeriodicBc();
	}
}

void CgnsBase::ReadAllCgnsZones( CgnsBase * cgnsBaseIn )
{
    cout << "** Reading CGNS Grid In Base " << this->baseId << "\n";
    cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
        cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
		CgnsZone * cgnsZone = this->cgnsZones[ iZone ];
        CgnsZone * cgnsZoneIn = cgnsBaseIn->cgnsZones[ iZone ];
		cgnsZone->ReadCgnsGrid( cgnsZoneIn );
	}
}

EndNameSpace