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
#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class FileIO;
class DataBook;

class ThermodynamicFunction
{
public:
    ThermodynamicFunction();
    ~ThermodynamicFunction();
public:
    int nTSpan, nPolyCoef;
    RealField2D polyCoef;
public:
    void Init( int nTSpan, int nPolyCoef );
    void ReadPolynomialCoefficient( FileIO * ioFile );
    void ReadPolynomialCoefficient( DataBook * dataBook );
    void WritePolynomialCoefficient( DataBook * dataBook );
};

class Thermodynamic
{
public:
    Thermodynamic();
    ~Thermodynamic();
public:
    int nSpecies;
    int nTSpan, nPolyCoef;
    HXVector< ThermodynamicFunction * > tfunction;
    RealField trange;
public:
    void Init( int nSpecies );
    void Read( FileIO * ioFile );
public:
    void Read( DataBook * dataBook );
    void Write( DataBook * dataBook );
public:
    RealField & GetPolyCoef( int is, int it );
    void GetTRangeId( const Real & dimt, int & itr );
//    Real GetTemperatureMin();
//    Real GetTemperatureMax();
};

EndNameSpace