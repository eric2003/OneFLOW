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

class SchmidtNumber;
class FileIO;
class DataBook;

class MolecularProperty
{
public:
    MolecularProperty();
    ~MolecularProperty();
public:
    int      nSpecies;
    RealField dim_mw;//dimensional molecular weight
    RealField mw; //molecular weight
    RealField dim_omw; 
    RealField omw;
    RealField ct; //characteristic temperature Of species;
    RealField mfrac; //reference mass fraction;
    StringField species_name;
    IntField ion_type;
    RealField cs; //collision cross section
    SchmidtNumber * schmidtNumber;
    Real dim_amw; //dimensional average molecular weight
    Real amw; //average molecular weight
public:
    void CalcProperty();
    void Init( int nSpecies );
    void Read( FileIO * ioFile );
    void Read( DataBook * dataBook );
    void Write( DataBook * dataBook );
};

EndNameSpace