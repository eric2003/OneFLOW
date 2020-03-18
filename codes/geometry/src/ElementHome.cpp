/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "HXPointer.h"

BeginNameSpace( ONEFLOW )

HXVector< UnitElement * > ElementHome::unitElement;
int ElementHome::numberOfUnitElement = 0;

ElementHome::ElementHome()
{
    ;
}

ElementHome::~ElementHome()
{
    ;
}

UnitElement * ElementHome::GetUnitElement( int elementType )
{
    return unitElement[ elementType ];
}

void ElementHome::Initialize()
{
    ElementHome::numberOfUnitElement = NofValidElementTypes;

    ONEFLOW::CreatePointer( unitElement, numberOfUnitElement );

    for ( UInt iUnitElement = 0; iUnitElement < numberOfUnitElement; ++ iUnitElement )
    {
        unitElement[ iUnitElement ]->Initialize( iUnitElement );
    }
}

void ElementHome::Free()
{
    ONEFLOW::DeletePointer( unitElement );
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
    UnitElement * unitElement = ElementHome::GetUnitElement( eType );
    int nodeNumber = unitElement->GetElementNodeNumbers( eType );
    return nodeNumber;
}

bool IsBasicVolumeElementType( int eType )
{
    return ElementHome::GetUnitElement( eType )->IsBasicVolumeElementType( eType );
}

EndNameSpace