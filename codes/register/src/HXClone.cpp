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

#include "HXClone.h"
#include "Fatal.h"
#include <map>
#include <memory>
#include <iostream>
#include <utility>


BeginNameSpace( ONEFLOW )

namespace
{
using CloneRegistry = std::map< std::string, std::unique_ptr< HXClone > >;

CloneRegistry & GetCloneRegistry()
{
    static CloneRegistry registry;
    return registry;
}
}

std::unique_ptr< HXClone > HXClone::SafeCloneUnique( const std::string & type )
{
    CloneRegistry & registry = GetCloneRegistry();
    CloneRegistry::iterator iter = registry.find( type );
    if ( iter == registry.end() )
    {
        Fatal( type + " class not found" );
        return nullptr;
    }

    return std::unique_ptr< HXClone >( iter->second->Clone() );
}


HXClone * HXClone::Register( const std::string & type, std::unique_ptr< HXClone > clone )
{
    //std::cout << "HXClone::Register : " << type << "\n";
    CloneRegistry & registry = GetCloneRegistry();
    CloneRegistry::iterator iter = registry.find( type );
    if ( iter != registry.end() ) return iter->second.get();

    HXClone * registeredClone = clone.get();
    registry.emplace( type, std::move( clone ) );
    return registeredClone;
}

HXClone * HXClone::Register( const std::string & type, HXClone * clone )
{
    return HXClone::Register( type, std::unique_ptr< HXClone >( clone ) );
}

EndNameSpace
