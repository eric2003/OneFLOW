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
#include "MolecularProperty.h"
#include "SchmidtNumber.h"
#include "FileIO.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

MolecularProperty::MolecularProperty()
{
    schmidtNumber = new SchmidtNumber();
}

MolecularProperty::~MolecularProperty()
{
    delete schmidtNumber;
}

void MolecularProperty::Init( int nSpecies )
{
    this->nSpecies = nSpecies;

    dim_mw.resize( nSpecies );
    dim_omw.resize( nSpecies );

    mw.resize( nSpecies );
    omw.resize( nSpecies );

    ct.resize( nSpecies );
    species_name.resize( nSpecies );

    ion_type.resize( nSpecies );
    mfrac.resize( nSpecies );
    cs.resize( nSpecies );

    schmidtNumber->Init( nSpecies );
}

void MolecularProperty::Read( FileIO * ioFile )
{
    //string word;
    string separator = " =\r\n#$,;\"'";
    ioFile->SetDefaultSeparator( separator );

    ioFile->SkipLines( 3 );
    //读各组元名称
    ioFile->ReadNextNonEmptyLine();
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        species_name[ iSpecies ] = ioFile->ReadNextWord();
    }

    ioFile->SkipLines( 3 );
    //read ion type
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        ion_type[ iSpecies ] = ioFile->ReadNextDigit< int >();
    }

    ioFile->SkipLines( 3 );
    //!读各组元分子量
    ioFile->ReadNextNonEmptyLine();
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dim_mw[ iSpecies ] = ioFile->ReadNextDigit< Real >();
    }

    ioFile->SkipLines( 3 );
    //read species mass fraction
    ioFile->ReadNextNonEmptyLine();
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        mfrac[ iSpecies ] = ioFile->ReadNextDigit< Real >();
    }

    ioFile->SkipLines( 3 );
    //read collision cross section
    ioFile->ReadNextNonEmptyLine();
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        cs[ iSpecies ] = ioFile->ReadNextDigit< Real >();
    }

    ioFile->SkipLines( 3 );
    //!读各组元特征温度
    ioFile->ReadNextNonEmptyLine();
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        ct[ iSpecies ] = ioFile->ReadNextDigit< Real >();
    }

    schmidtNumber->CalcSchmidtNumber( ion_type );
}

void MolecularProperty::Read( DataBook * dataBook )
{
    schmidtNumber->Read( dataBook );
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dataBook->ReadString( species_name[ iSpecies ] );
    }

    HXRead( dataBook, ion_type );
    HXRead( dataBook, mfrac );

    HXRead( dataBook, cs );
    HXRead( dataBook, ct );

    HXRead( dataBook, dim_mw );
    HXRead( dataBook, mw );

    HXRead( dataBook, dim_omw );
    HXRead( dataBook, omw );
}

void MolecularProperty::Write( DataBook * dataBook )
{
    schmidtNumber->Write( dataBook );
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dataBook->AppendString( species_name[ iSpecies ] );
    }

    HXAppend( dataBook, ion_type );
    HXAppend( dataBook, mfrac );

    HXAppend( dataBook, cs );
    HXAppend( dataBook, ct );

    HXAppend( dataBook, dim_mw );
    HXAppend( dataBook, mw );

    HXAppend( dataBook, dim_omw );
    HXAppend( dataBook, omw );
}

void MolecularProperty::CalcProperty()
{
    //dimensional molecular weight
    Real coef = 1.0e-3;
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dim_mw[ iSpecies ] *= coef;
    }

    //dimensional molecular weight ( reciprocal )
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dim_omw[ iSpecies ] = 1.0 / dim_mw[ iSpecies ];
    }

    //dimensional average molecular weight( reciprocal )
    Real dim_oamw = 0.0;
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        dim_oamw += mfrac[ iSpecies ] * dim_omw[ iSpecies ];
    }

    //dimensional average molecular weight
    dim_amw = 1.0 / dim_oamw;

    //average molecular weight
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        mw[ iSpecies ] = dim_mw[ iSpecies ] * dim_oamw;
    }

    //average molecular weight( reciprocal )
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        omw[ iSpecies ] = 1.0 / mw[ iSpecies ];
    }

    amw = 1.0;
}


EndNameSpace