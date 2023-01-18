/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include <iostream>


BeginNameSpace( ONEFLOW )

std::map< std::string, HXClone * > * HXClone::classMap = 0;

HXClone * HXClone::SafeClone( const std::string & type )
{
    std::map < std::string, HXClone * >::iterator iter = HXClone::classMap->find( type );
    if ( iter == HXClone::classMap->end() )
    {
        std::cout << type << " class not found" << std::endl;
        exit( 0 );
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
