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
#include "Thermodynamic.h"
#include "FileIO.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

ThermodynamicFunction::ThermodynamicFunction()
{
    ;
}

ThermodynamicFunction::~ThermodynamicFunction()
{
    ;
}

void ThermodynamicFunction::Init( int nTSpan, int nPolyCoef )
{
    this->nTSpan = nTSpan;
    this->nPolyCoef = nPolyCoef;
    AllocateVector( polyCoef, nTSpan, nPolyCoef );
}

void ThermodynamicFunction::ReadPolynomialCoefficient( FileIO * ioFile )
{
    string word;
    string separator = " =\r\n#$,;\"'";

    ioFile->SetDefaultSeparator( separator );

    for ( int iInterval = 0; iInterval < nTSpan; ++ iInterval )
    {
        ioFile->ReadNextNonEmptyLine();

        RealField & coef = polyCoef[ iInterval ];

        for ( int icoef = 0; icoef < nPolyCoef; ++ icoef )
        {
            coef[ icoef ] = ioFile->ReadNextDigit< Real >();
        }
    }
}

void ThermodynamicFunction::ReadPolynomialCoefficient( DataBook * dataBook )
{
    for ( int iInterval = 0; iInterval < nTSpan; ++ iInterval )
    {
        RealField & coef = polyCoef[ iInterval ];

        HXRead( dataBook, coef );
    }
}

void ThermodynamicFunction::WritePolynomialCoefficient( DataBook * dataBook )
{
    for ( int iInterval = 0; iInterval < nTSpan; ++ iInterval )
    {
        RealField & coef = polyCoef[ iInterval ];

        HXAppend( dataBook, coef );
    }
}


Thermodynamic::Thermodynamic()
{
    ;
}

Thermodynamic::~Thermodynamic()
{
    size_t nSize = tfunction.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete tfunction[ i ];
    }
}

void Thermodynamic::Init( int nSpecies )
{
    this->nSpecies = nSpecies;
    tfunction.resize( nSpecies );
    for ( int i = 0; i < nSpecies; ++ i )
    {
        tfunction[ i ] = new ThermodynamicFunction();
    }
}

void Thermodynamic::Read( FileIO * ioFile )
{
    string word;
    string separator = " =\r\n#$,;\"'";

    ioFile->SetDefaultSeparator( separator );

    ioFile->SkipLines( 3 );
    ioFile->ReadNextNonEmptyLine();

    //read number of polynomial coefficient.( a1-a7, or a1-a6 )
    nTSpan = ioFile->ReadNextDigit< int >();
    nPolyCoef = ioFile->ReadNextDigit< int >();

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        tfunction[ iSpecies ]->Init( nTSpan, nPolyCoef );
    }

    trange.resize( nTSpan + 1 );

    ioFile->SkipLines( 4 );
    ioFile->ReadNextNonEmptyLine();

    for ( int m = 0; m < nTSpan + 1; ++ m )
    {
        trange[ m ] = ioFile->ReadNextDigit< Real >();
    }

    ioFile->SkipLines( 2 );

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        ioFile->SkipLines( 3 );
        tfunction[ iSpecies ]->ReadPolynomialCoefficient( ioFile );
    }
}

void Thermodynamic::Read( DataBook * dataBook )
{
    HXRead( dataBook, nTSpan );
    HXRead( dataBook, nPolyCoef );

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        tfunction[ iSpecies ]->Init( nTSpan, nPolyCoef );
    }

    trange.resize( nTSpan + 1 );

    HXRead( dataBook, trange );
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        tfunction[ iSpecies ]->ReadPolynomialCoefficient( dataBook );
    }
}

void Thermodynamic::Write( DataBook * dataBook )
{
    HXAppend( dataBook, nTSpan );
    HXAppend( dataBook, nPolyCoef );
    HXAppend( dataBook, trange );

    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        tfunction[ iSpecies ]->WritePolynomialCoefficient( dataBook );
    }
}

RealField & Thermodynamic::GetPolyCoef( int is, int it )
{
    return tfunction[ is ]->polyCoef[ it ];
}

void Thermodynamic::GetTRangeId( const Real & dimt, int & itr )
{
    itr = 0;
    if ( trange[ 0 ] >= dimt )
    {
        itr = 0;
    }
    else if ( trange[ nTSpan ] <= dimt )
    {
        itr = nTSpan - 1;
    }
    else
    {
        for ( int it = 0; it < nTSpan; ++ it )
        {
            if ( ( trange[ it ] < dimt ) &&
                 ( dimt <= trange[ it + 1 ] ) )
            {
                itr = it;
                break;
            }
        }
    }
}

EndNameSpace