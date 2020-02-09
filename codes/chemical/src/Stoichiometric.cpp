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
#include "Stoichiometric.h"
#include "FileIO.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

Stoichiometric::Stoichiometric()
{
    ;
}

Stoichiometric::~Stoichiometric()
{
    ;
}

void Stoichiometric::Init( int nSpecies, int nReaction )
{
    this->nSpecies  = nSpecies;
    this->nReaction = nReaction;

    AllocateVector( mf, nReaction, nSpecies );
    AllocateVector( mb, nReaction, nSpecies );
    AllocateVector( mt, nReaction, nSpecies );
}

void Stoichiometric::Read( FileIO * ioFile )
{
    string word;
    string separator = " =\r\n#$,;\"'";

    ioFile->SetDefaultSeparator( separator );

    if ( nReaction <= 0 ) return;

    ioFile->SkipLines( 4 );
    //forward reaction stoichiometric coefficient matrix
    for ( int iReaction = 0; iReaction < nReaction; ++ iReaction )
    {
        ioFile->ReadNextNonEmptyLine();
        //¶Áirtmp
        word = ioFile->ReadNextWord();

        for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
        {
            mf[ iReaction ][ iSpecies ] = ioFile->ReadNextDigit< int >();
        }
    }

    ioFile->SkipLines( 4 );
    //bakward reaction stoichiometric coefficient matrix
    for ( int iReaction = 0; iReaction < nReaction; ++ iReaction )
    {
        ioFile->ReadNextNonEmptyLine();
        //irtmp
        word = ioFile->ReadNextWord();

        for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
        {
            mb[ iReaction ][ iSpecies ] = ioFile->ReadNextDigit< int >();
        }
    }

    ioFile->SkipLines( 4 );
    //third body reaction stoichiometric coefficient matrix
    for ( int iReaction = 0; iReaction < nReaction; ++ iReaction )
    {
        ioFile->ReadNextNonEmptyLine();
        //irtmp
        word = ioFile->ReadNextWord();

        for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
        {
            mt[ iReaction ][ iSpecies ] = ioFile->ReadNextDigit< Real >();
        }
    }
}

void Stoichiometric::Read( DataBook * dataBook )
{
    HXRead( dataBook, mf );
    HXRead( dataBook, mb );
    HXRead( dataBook, mt );
}

void Stoichiometric::Write( DataBook * dataBook )
{
    HXAppend( dataBook, mf );
    HXAppend( dataBook, mb );
    HXAppend( dataBook, mt );
}

EndNameSpace