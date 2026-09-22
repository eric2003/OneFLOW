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

#include "UsdPara.h"
#include <iostream>


BeginNameSpace( ONEFLOW )


UsdPara::UsdPara()
{
    this->flag = false;
}

UsdPara::~UsdPara()
{
    ;
}

void UsdPara::Init(
    const UsdFieldNames & fieldNames,
    int nEqu )
{
    if ( this->flag ) return;

    this->flag = true;

    this->flow.push_back( fieldNames.q );
    this->flow.push_back( fieldNames.q1 );
    this->flow.push_back( fieldNames.q2 );

    this->residual.push_back( fieldNames.res );
    this->residual.push_back( fieldNames.res1 );
    this->residual.push_back( fieldNames.res2 );

    this->dq.push_back( fieldNames.dq );

    this->nEqu = nEqu;
}

EndNameSpace
