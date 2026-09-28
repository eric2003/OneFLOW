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
#include <iostream>


BeginNameSpace( ONEFLOW )

std::map< std::string, HXClone * > * HXClone::classMap = 0;

HXClone * HXClone::SafeClone( const std::string & type )
{
    // FIX: classMap may be null if nothing has been Register()'d yet.
    // Fatal(...) throws std::runtime_error, so callers of SafeClone must
    // be prepared to handle that exception (or let it propagate) rather
    // than expecting a null return - the `return nullptr;` lines below
    // are unreachable and exist only to satisfy the compiler.
    if ( ! HXClone::classMap )
    {
        Fatal( type + " class not found" );
        return nullptr;
    }

    auto iter = HXClone::classMap->find( type );
    if ( iter == HXClone::classMap->end() )
    {
        Fatal( type + " class not found" );
        return nullptr;
    }

    return iter->second->Clone();
}

HXClone * HXClone::Register( const std::string & type, HXClone * clone )
{
    if ( ! HXClone::classMap )
    {
        HXClone::classMap = new std::map < std::string, HXClone * >();
    }

    //std::cout << "HXClone::Register : " << type << "\n";

    std::map < std::string, HXClone * >::iterator iter = HXClone::classMap->find( type );
    if ( iter == HXClone::classMap->end() )
    {
        ( * HXClone::classMap )[ type ] = clone;
        return clone;
    }
    else
    {
        delete clone;
        return iter->second;
    }
}

EndNameSpace
