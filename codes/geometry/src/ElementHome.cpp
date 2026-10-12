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

BeginNameSpace( ONEFLOW )

HXVector< std::unique_ptr< UnitElement > > ElementHome::unitElement;
int ElementHome::numberOfUnitElement = 0;

UnitElement & ElementHome::GetUnitElement( int elementType )
{
    return *unitElement[ elementType ];
}

void ElementHome::Initialize()
{
    ElementHome::numberOfUnitElement = NofValidElementTypes;

    unitElement.clear();
    unitElement.reserve( numberOfUnitElement );

    for ( HXSize_t iUnitElement = 0; iUnitElement < numberOfUnitElement; ++ iUnitElement )
    {
        auto unit = std::make_unique< UnitElement >();
        unit->Initialize( iUnitElement );
        unitElement.push_back( std::move( unit ) );
    }
}

void ElementHome::Free()
{
    unitElement.clear();
    ElementHome::numberOfUnitElement = 0;
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
    UnitElement & unitElement = ElementHome::GetUnitElement( eType );
    int nodeNumber = unitElement.GetElementNodeNumbers( eType );
    return nodeNumber;
}

int GetElementDimension( int eType )
{
    if ( eType < 0 || eType >= NofValidElementTypes ) return -1;

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
    return ElementHome::GetUnitElement( eType ).IsBasicVolumeElementType( eType );
}

EndNameSpace
