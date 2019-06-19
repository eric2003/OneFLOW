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

CgnsFamilyBc::CgnsFamilyBc()
{
	int_flag = false;
	Init();
}

CgnsFamilyBc::~CgnsFamilyBc()
{
	Free();
}

void CgnsFamilyBc::Init()
{
	if ( int_flag ) return;
	int_flag = true;
	bcMap = new map< string, int >;
}

void CgnsFamilyBc::Free()
{
	delete bcMap;
}

void CgnsFamilyBc::Register( const string & regionName, int bcType )
{
	map< string, int >::iterator iter = bcMap->find( regionName );
	if ( iter == bcMap->end() )
	{
		( * CgnsFamilyBc::bcMap )[ regionName ] = bcType;
	}
}

void CgnsFamilyBc::Unregister( const string & regionName )
{
	bcMap->erase( regionName );
}

int CgnsFamilyBc::GetBcType( const string & regionName )
{
	map< string, int >::iterator iter = bcMap->find( regionName );
	if ( iter == bcMap->end() )
	{
		return -1;
	}

	return iter->second;
}


EndNameSpace