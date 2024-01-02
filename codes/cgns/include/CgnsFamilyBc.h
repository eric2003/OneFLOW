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
#include "HXCgns.h"
#include <string>
#include <map>

BeginNameSpace( ONEFLOW )


#ifdef ENABLE_CGNS

class CgnsBase;

class CgnsFamilyBc
{
public:
    CgnsFamilyBc( CgnsBase * cgnsBase );
    ~CgnsFamilyBc();
public:
    std::map< std::string, int > * bcMap;
    CgnsBase * cgnsBase;
public:
    void Init();
    void Free();
    void Register( const std::string & regionName, int bcType );
    void Unregister( const std::string & regionName );
    int GetBcType( const std::string & regionName );
public:
    void SetFamilyBc( BCType_t & bcType, const std::string & bcRegionName );
    BCType_t GetFamilyBcType( const std::string & bcFamilyName );
    void ReadFamilySpecifiedBc();
};

#endif

EndNameSpace
