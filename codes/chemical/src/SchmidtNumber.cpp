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
#include "SchmidtNumber.h"
#include "HXMath.h"
#include "NsCom.h"
#include "DataBook.h"
#include "DataBaseIO.h"

BeginNameSpace( ONEFLOW )

SchmidtNumber::SchmidtNumber()
{
    ;
}

SchmidtNumber::~SchmidtNumber()
{
    ;
}

void SchmidtNumber::Init( int nSpecies )
{
    this->nSpecies = nSpecies;

    lamSchmidt.resize( nSpecies );
    turbSchmidt.resize( nSpecies );
    olamSchmidt.resize( nSpecies );
    oturbSchmidt.resize( nSpecies );
}

void SchmidtNumber::CalcSchmidtNumber( IntField & ionType )
{
    for ( int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies )
    {
        if ( ionType[ iSpecies ] != 0 && ABS( ionType[ iSpecies ] ) < 100 )
        {
            //if ( ionType[ iSpecies ] > 0 )
            //{
            lamSchmidt [ iSpecies ] = half * nscom.schmidtl;
            turbSchmidt[ iSpecies ] = half * nscom.schmidtt;
            //}
            //else
            //{
            //    lamSchmidt [ iSpecies ] = half * nscom.schmidtl;
            //    turbSchmidt[ iSpecies ] = half * nscom.schmidtt;
            //}
        }
        else
        {
            lamSchmidt [ iSpecies ] = 1.0 * nscom.schmidtl;
            turbSchmidt[ iSpecies ] = 1.0 * nscom.schmidtt;
        }

        olamSchmidt [ iSpecies ] = 1.0 / lamSchmidt [ iSpecies ];
        oturbSchmidt[ iSpecies ] = 1.0 / turbSchmidt[ iSpecies ];
    }
}

void SchmidtNumber::Read( DataBook * dataBook )
{
    HXRead( dataBook, lamSchmidt );
    HXRead( dataBook, turbSchmidt );

    HXRead( dataBook, olamSchmidt );
    HXRead( dataBook, oturbSchmidt );
}

void SchmidtNumber::Write( DataBook * dataBook )
{
    HXAppend( dataBook, lamSchmidt );
    HXAppend( dataBook, turbSchmidt );

    HXAppend( dataBook, olamSchmidt );
    HXAppend( dataBook, oturbSchmidt );
}


EndNameSpace