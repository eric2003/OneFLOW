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

#include "ScalarCgns.h"

BeginNameSpace( ONEFLOW )

SectionMarker::SectionMarker()
{
    ;
}

SectionMarker::~SectionMarker()
{
    ;
}

SectionManager::SectionManager()
{
    ;
}

SectionManager::~SectionManager() = default;

int SectionManager::GetNSections()
{
    return this->data.size();
}

void SectionManager::Alloc( int nType )
{
    this->nType = nType;
    this->data.resize( nType );
    for ( int i = 0; i < nType; ++ i )
    {
        this->data[ i ] = std::make_unique< SectionMarker >();
    }
}

int SectionManager::CalcTotalElem()
{
    int nType = this->data.size();
    int nElements = 0;
    for ( int i = 0; i < nType; ++ i )
    {
        SectionMarker * sectionMarker = this->data[ i ].get();
        nElements += sectionMarker->nElements;
    }
    return nElements;
}

EndNameSpace
