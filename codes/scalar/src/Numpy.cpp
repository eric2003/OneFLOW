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

#include "Numpy.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

Numpy::Numpy()
{
}

Numpy::~Numpy()
{
}

void Numpy::Ones( vector< double > & var )
{
    int nSize = var.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        var[ i ] = 1.0;
    }
}

void Numpy::Set( vector< double > & var, int st, int ed, double v )
{
    for ( int i = st; i < ed; ++ i )
    {
        var[ i ] = v;
    }
}

void Numpy::Linspace( vector< double > & var, double st, double ed )
{
    int nSize = var.size();
    double ds = ed - st;
    double dx = ds / ( nSize - 1.0 );
    for ( int i = 0; i < nSize; ++ i )
    {
        var[ i ] = st + i * dx;
    }
}

void Numpy::Plot( vector< double > & x, vector< double > & f )
{
    fstream file;
    file.open( "plot.plt", ios_base::out );
    int nSize = x.size();
    file << nSize << " ";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << x[ i ] << " " << f[ i ] << " ";
    }

    file.close();
    file.clear();
}


EndNameSpace