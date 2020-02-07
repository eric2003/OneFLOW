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
#include "CgnsGlobal.h"
#include "CgnsMultiBase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsGlobal cgns_global;

CgnsGlobal::CgnsGlobal()
{
    ;
}

CgnsGlobal::~CgnsGlobal()
{
    ;
}

CgnsZone * CgnsGlobal::GetCgnsZone( int zId )
{
    return cgnsbases->GetCgnsZone( zId - 1 );
}

CgnsZone * CgnsGlobal::GetCgnsZone( const string & zoneName )
{
    CgnsBase * cgnsBase = cgnsbases->baseVector[ 0 ];
    return cgnsBase->GetCgnsZone( zoneName );
}

CgnsZone * GetCgnsZone( const string & zoneName )
{
    return cgns_global.GetCgnsZone( zoneName );
}

#endif
EndNameSpace