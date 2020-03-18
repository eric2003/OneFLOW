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

#include "UsdPara.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UsdPara::UsdPara()
{
    this->flag = false;
}

UsdPara::~UsdPara()
{
    ;
}

void UsdPara::Init( StringField & fieldNameList, int nEqu )
{
    if ( this->flag ) return;
    this->flag = true;

    this->flow.push_back( fieldNameList[ 0 ] );
    this->flow.push_back( fieldNameList[ 1 ] );
    this->flow.push_back( fieldNameList[ 2 ] );

    this->residual.push_back( fieldNameList[ 3 ] );
    this->residual.push_back( fieldNameList[ 4 ] );
    this->residual.push_back( fieldNameList[ 5 ] );

    this->dq.push_back( fieldNameList[ 6 ] );

    this->nEqu = nEqu;
}

EndNameSpace