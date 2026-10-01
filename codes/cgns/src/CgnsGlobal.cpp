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
#include "CgnsGlobal.h"
#include "CgnsZbase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include <stdexcept>

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsGlobal cgns_global;

CgnsGlobal::CgnsGlobal()
    : cgnsbases( nullptr )
{
}

void CgnsGlobal::Bind( CgnsZbase * cgnsBases )
{
    cgnsbases = cgnsBases;
}

void CgnsGlobal::ClearIfBoundTo( const CgnsZbase * cgnsBases )
{
    if ( cgnsbases == cgnsBases )
    {
        cgnsbases = nullptr;
    }
}

bool CgnsGlobal::IsBoundTo( const CgnsZbase * cgnsBases ) const
{
    return cgnsbases == cgnsBases;
}

CgnsGlobal::~CgnsGlobal()
{
    ;
}

CgnsZone * CgnsGlobal::GetCgnsZoneByName( const std::string & zoneName )
{
    if ( cgnsbases == nullptr )
    {
        throw std::logic_error( "CGNS zone lookup requested without an active CGNS base" );
    }
    CgnsBase * cgnsBase = cgnsbases->baseVector[ 0 ].get();
    return cgnsBase->GetCgnsZoneByName( zoneName );
}

CgnsZone * GetCgnsZoneByName( const std::string & zoneName )
{
    return cgns_global.GetCgnsZoneByName( zoneName );
}

#endif
EndNameSpace
