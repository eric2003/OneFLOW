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

#include "Scalar.h"
#include "Numpy.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

Scalar::Scalar()
{
}

Scalar::~Scalar()
{
}

void Scalar::Run()
{
    int nx = 41;

    double len = 2.0;
    double dx = len / ( nx - 1.0 );
    int    nt = 25;
    double dt = 0.025;
    double c  = 1;

    vector< double > u( nx );
    Numpy::Ones( u );

    int st = 0.5 / dx;
    int ed = 1 / dx + 1;

    Numpy::Set( u, st, ed, 2 );

    vector< double > x( nx );
    Numpy::Linspace( x, 0, len );
    Numpy::Plot( x, u );
}


EndNameSpace