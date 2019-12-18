/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

class Stoichiometric
{
public:
    Stoichiometric();
    ~Stoichiometric();
public:
    int nSpecies, nReaction;
    //forward and bakward reaction stoichiometric coefficient matrix
    LinkField mf, mb;
    //third body reaction stoichiometric coefficient matrix;
    RealField2D mt;
public:
    void Init( int nSpecies, int nReaction );
    void Read( FileIO * ioFile );
    void Read ( DataBook * dataBook );
    void Write( DataBook * dataBook );
};

EndNameSpace