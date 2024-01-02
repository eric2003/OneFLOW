/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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


#pragma once
#include "HXVector.h"
#include <string>

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS
#include "cgnslib.h"

typedef cgsize_t CgInt;
typedef HXVector< CgInt > CgIntField;
typedef HXVector< CgIntField > CgLinkField;

class CgnsTraits
{
public:
    typedef char char33[ 33 ];
};

std::string GetCgnsPointSetName         ( int cgnsPointSetType );
std::string GetCgnsBcName( int cgnsBcType );
std::string GetCgnsGridLocationName     ( int cgnsGridLocation );
std::string GetCgnsZoneTypeName         ( int cgnsZoneType );

#endif
EndNameSpace
