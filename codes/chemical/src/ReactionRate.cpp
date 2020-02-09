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
#include "ReactionRate.h"
#include "FileIO.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

ReactionRate::ReactionRate()
{
    ;
}

ReactionRate::~ReactionRate()
{
    ;
}

void ReactionRate::Init( int nReaction )
{
    this->nReaction = nReaction;

    f1.resize( nReaction );
    f2.resize( nReaction );
    f3.resize( nReaction );

    b1.resize( nReaction );
    b2.resize( nReaction );
    b3.resize( nReaction );
}

void ReactionRate::Read( FileIO * ioFile )
{
    string separator = " =\r\n#$,;\"'";

    ioFile->SetDefaultSeparator( separator );

    if ( nReaction <= 0 ) return;

    ioFile->SkipLines( 3 );
    for ( int iReaction = 0; iReaction < nReaction; ++ iReaction )
    {
        ioFile->ReadNextNonEmptyLine();
        //¶Áirtmp
        string word = ioFile->ReadNextWord();

        f1[ iReaction ] = ioFile->ReadNextDigit< Real >();
        f2[ iReaction ] = ioFile->ReadNextDigit< Real >();
        f3[ iReaction ] = ioFile->ReadNextDigit< Real >();

        b1[ iReaction ] = ioFile->ReadNextDigit< Real >();
        b2[ iReaction ] = ioFile->ReadNextDigit< Real >();
        b3[ iReaction ] = ioFile->ReadNextDigit< Real >();
    }
}

void ReactionRate::Read( DataBook * dataBook )
{
    HXRead( dataBook, f1 );
    HXRead( dataBook, f2 );
    HXRead( dataBook, f3 );

    HXRead( dataBook, b1 );
    HXRead( dataBook, b2 );
    HXRead( dataBook, b3 );
}

void ReactionRate::Write( DataBook * dataBook )
{
    HXAppend( dataBook, f1 );
    HXAppend( dataBook, f2 );
    HXAppend( dataBook, f3 );

    HXAppend( dataBook, b1 );
    HXAppend( dataBook, b2 );
    HXAppend( dataBook, b3 );
}

EndNameSpace