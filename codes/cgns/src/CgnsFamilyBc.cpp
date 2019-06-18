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

#include "CgnsFamilyBc.h"

BeginNameSpace( ONEFLOW )

bool FamilyBc::int_flag = false;
map< string, int > * FamilyBc::bcMap = 0;

FamilyBc::FamilyBc()
{
	;
}

FamilyBc::~FamilyBc()
{
	;
}

void FamilyBc::Init()
{
	if ( FamilyBc::int_flag ) return;
	FamilyBc::int_flag = true;
	FamilyBc::bcMap = new map< string, int >;
}

void FamilyBc::Free()
{
	delete FamilyBc::bcMap;
}

void FamilyBc::Register( const string & regionName, int bcType )
{
	map< string, int >::iterator iter = FamilyBc::bcMap->find( regionName );
	if ( iter == FamilyBc::bcMap->end() )
	{
		( * FamilyBc::bcMap )[ regionName ] = bcType;
	}
}

void FamilyBc::Unregister( const string & regionName )
{
	FamilyBc::bcMap->erase( regionName );
}

int FamilyBc::GetBcType( const string & regionName )
{
	map< string, int >::iterator iter = FamilyBc::bcMap->find( regionName );
	if ( iter == FamilyBc::bcMap->end() )
	{
		return -1;
	}

	return iter->second;
}


EndNameSpace