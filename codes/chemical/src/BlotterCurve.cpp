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
#include "BlotterCurve.h"
#include "FileIO.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

BlotterCurve::BlotterCurve()
{
    ;
}

BlotterCurve::~BlotterCurve()
{
    ;
}

void BlotterCurve::Init( int nSpecies )
{
    this->nSpecies = nSpecies;

    a.resize( nSpecies );
    b.resize( nSpecies );
    c.resize( nSpecies );
    d.resize( nSpecies );
    e.resize( nSpecies );
}

void BlotterCurve::Read( FileIO * ioFile )
{
    string word;
    string separator = " =\r\n#$,;\"'";

    ioFile->SetDefaultSeparator( separator );

    ioFile->SkipLines( 2 );

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        ioFile->SkipLines( 3 );
        ioFile->ReadNextNonEmptyLine();

        a[ iSpecies ] = ioFile->ReadNextDigit< Real >();
        b[ iSpecies ] = ioFile->ReadNextDigit< Real >();
        c[ iSpecies ] = ioFile->ReadNextDigit< Real >();
        d[ iSpecies ] = ioFile->ReadNextDigit< Real >();
        e[ iSpecies ] = ioFile->ReadNextDigit< Real >();
    }
}

void BlotterCurve::Read( DataBook * dataBook )
{
    HXRead( dataBook, a );
    HXRead( dataBook, b );
    HXRead( dataBook, c );
    HXRead( dataBook, d );
    HXRead( dataBook, e );
}

void BlotterCurve::Write( DataBook * dataBook )
{
    HXAppend( dataBook, a );
    HXAppend( dataBook, b );
    HXAppend( dataBook, c );
    HXAppend( dataBook, d );
    HXAppend( dataBook, e );
}

EndNameSpace