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
#pragma once
#include <memory>
#include "HXDefine.h"
#include "DataBase.h"
#include "DataPointer.h"
#include <string>

BeginNameSpace( ONEFLOW )

template< typename T >
void CreateMRField( T * storage, int nEqu, int nSize, const std::string & fieldName )
{
    auto mrField = std::make_unique<MRField>( nEqu, nSize );
    auto wrap = std::make_unique<DataPointer<MRField>>( std::move( mrField ) );
    ONEFLOW::CreateFieldPointer( storage, std::move( wrap ), fieldName );
}

void ZeroField( MRField * field, int nEqu, int nSize );

EndNameSpace
