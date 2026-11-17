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


#include "DataBaseIO.h"
#include <vector>
#include "DataBook.h"

BeginNameSpace( ONEFLOW )

void HXRead( DataBook * dataBook, std::string & cs )
{
    int nLength = 0;
    ONEFLOW::HXRead( dataBook, nLength );

    std::vector<char> data( static_cast<std::size_t>( nLength ) + 1, '\0' );
    dataBook->Read( data.data(), nLength + 1 );

    cs = data.data();
}

void HXWrite( DataBook * dataBook, const std::string & cs )
{
    int nLength = static_cast<int>( cs.length() );
    ONEFLOW::HXWrite( dataBook, nLength );

    std::vector<char> data( static_cast<std::size_t>( nLength ) + 1, '\0' );
    cs.copy( data.data(), nLength );

    dataBook->Write( data.data(), nLength + 1 );
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
