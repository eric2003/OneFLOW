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


#include "DataBaseIO.h"
#include "DataBook.h"

BeginNameSpace( ONEFLOW )

void HXRead( DataBook * dataBook, string & cs )
{
    int nLength = 0;
    ONEFLOW::HXRead( dataBook, nLength );

    char * data = new char[ nLength + 1 ];
    dataBook->Read( data, nLength + 1 );

    cs = data;

    delete[] data;
}

void HXWrite( DataBook * dataBook, string & cs )
{
    int nLength = cs.length();
    ONEFLOW::HXWrite( dataBook, nLength );

    char * data = new char[ nLength + 1 ];

    cs.copy( data, nLength );
    data[ nLength ] = '\0';

    dataBook->Write( data, nLength + 1 );

    delete[] data;
}

void HXRead( DataBook * dataBook, MRField * field )
{
    int nEqu = field->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXRead( dataBook, ( * field )[ iEqu ] );
    }
}

void HXWrite( DataBook * dataBook, MRField * field )
{
    int nEqu = field->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXWrite( dataBook, ( * field )[ iEqu ] );
    }
}

EndNameSpace
