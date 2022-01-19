/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

#include "ScalarDataIO.h"
#include "DataBase.h"
#include "ActionState.h"


BeginNameSpace( ONEFLOW )

void HXWriteField( DataBook * dataBook, MRField * field2D, std::vector< int > & idMap )
{
    int nEqu = field2D->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXWriteField( dataBook, ( * field2D )[ iEqu ], idMap );
    }
}

void HXWriteField( DataBook * dataBook, RealField & field, std::vector< int > & idMap )
{
    int nElem = idMap.size();
    if ( nElem <= 0 ) return;
    RealField swapField( nElem );
    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int id = idMap[ iElem ];
        swapField[ iElem ] = field[ id ];
    }
    HXWrite( dataBook, swapField );
}

void HXReadField( DataBook * dataBook, MRField * field2D, std::vector< int > & idMap )
{
    int nEqu = field2D->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXReadField( dataBook, ( * field2D )[ iEqu ], idMap );
    }
}

void HXReadField( DataBook * dataBook, RealField & field, std::vector< int > &  idMap )
{
    int nElem = idMap.size();
    if ( nElem <= 0 ) return;
    RealField swapField( nElem );
    HXRead( dataBook, swapField );

    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int id = idMap[ iElem ];
        field[ id ] = swapField[ iElem ];
    }
}

EndNameSpace
