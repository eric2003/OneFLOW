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

#include "ElementHome.h"
#include "HXCgns.h"
#include <utility>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

HXVector< std::unique_ptr< const UnitElement > > ElementHome::unitElements;

int ElementHome::GetNumberOfUnitElements()
{
    return static_cast< int >( unitElements.size() );
}

bool ElementHome::HasUnitElement( int elementType )
{
    return elementType >= 0 &&
        static_cast< std::size_t >( elementType ) < unitElements.size() &&
        static_cast< bool >( unitElements[ static_cast< std::size_t >( elementType ) ] );
}

const UnitElement & ElementHome::GetUnitElementDefinition( int elementType )
{
    if ( ! HasUnitElement( elementType ) )
    {
        throw std::out_of_range( "ElementHome::GetUnitElementDefinition: element type is out of range" );
    }

    return *unitElements[ static_cast< std::size_t >( elementType ) ];
}

void ElementHome::Initialize()
{
    // Keep enum values as registry indices, including unsupported gaps.
    // A null slot means that the element type has no UnitElement definition.
    HXVector< std::unique_ptr< const UnitElement > > initializedElements;
    initializedElements.resize( NofValidElementTypes );

    for ( HXSize_t iUnitElement = 0; iUnitElement < NofValidElementTypes; ++ iUnitElement )
    {
        auto unit = std::make_unique< UnitElement >();
        if ( unit->GetElementNodeNumbers( iUnitElement ) <= 0 )
        {
            continue;
        }

        unit->Initialize( iUnitElement );
        initializedElements[ iUnitElement ] = std::move( unit );
    }

    unitElements.swap( initializedElements );
}

void ElementHome::Free()
{
    unitElements.clear();
}

class ElementHomeInit
{
public:
    ElementHomeInit()
    {
        ElementHome::Initialize();
    }
    ~ElementHomeInit()
    {
        ElementHome::Free();
    }
};
ElementHomeInit elementHomeInit;

int GetElementNodeNumbers( int eType )
{
    if ( eType < 0 || eType >= ElementHome::GetNumberOfUnitElements() )
    {
        throw std::out_of_range( "GetElementNodeNumbers: element type is out of range" );
    }

    if ( ! ElementHome::HasUnitElement( eType ) ) return 0;

    return ElementHome::GetUnitElementDefinition( eType ).GetElementNodeNumbers( eType );
}

int GetElementDimension( int eType )
{
    if ( ! ElementHome::HasUnitElement( eType ) ) return -1;

    switch ( eType )
    {
    case NODE:
        return 0;
    case BAR_2:
    case BAR_3:
        return 1;
    case TRI_3:
    case TRI_6:
    case QUAD_4:
    case QUAD_8:
    case QUAD_9:
        return 2;
    case TETRA_4:
    case TETRA_10:
    case PYRA_5:
    case PYRA_14:
    case PENTA_6:
    case PENTA_15:
    case PENTA_18:
    case HEXA_8:
    case HEXA_20:
    case HEXA_27:
        return 3;
    default:
        return -1;
    }
}

bool IsBasicVolumeElementType( int eType )
{
    if ( ! ElementHome::HasUnitElement( eType ) ) return false;

    return ElementHome::GetUnitElementDefinition( eType ).IsBasicVolumeElementType( eType );
}

EndNameSpace
