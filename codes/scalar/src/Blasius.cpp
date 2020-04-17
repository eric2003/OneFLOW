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

#include "Blasius.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

Blasius::Blasius()
{
}

Blasius::~Blasius()
{
}

void Blasius::Run()
{
    double deta = 0.01;
    this->rk4( deta );
}

void Blasius::BlasiusFun( vector<double> &y, vector<double> &k )
{
    k[ 0 ] = y[ 1 ];
    k[ 1 ] = y[ 2 ];
    k[ 2 ] = - 0.5 * y[ 0 ] * y[ 2 ];
}

void Blasius::rk4( double deta )
{
    int nvec = 3;
    vector< double > y1( nvec, 0 ), y2( nvec, 0 ), y3( nvec, 0 ), y4( nvec, 0 );
    vector< double > k1( nvec, 0 ), k2( nvec, 0 ), k3( nvec, 0 ), k4( nvec, 0 );
    double eta = 0.0;
    double etaN = 10.0;
    int n = int( etaN / deta );
    double h = deta;

    y1[ 0 ] = 0.0;
    y1[ 1 ] = 0.0;
    y1[ 2 ] = 1.0;

    double coef = 1.0 / 6.0 * h;

    for ( int i = 0; i < n; ++ i )
    {
        BlasiusFun( y1, k1 );
        for ( int m = 0; m < nvec; ++ m )
        {
            y2[ m ] = y1[ m ] + 0.5 * h * k1[ m ];
        }

        BlasiusFun( y2, k2 );

        for ( int m = 0; m < nvec; ++ m )
        {
            y3[ m ] = y1[ m ] + 0.5 * h * k2[ m ];
        }

        BlasiusFun( y3, k3 );

        for ( int m = 0; m < nvec; ++ m )
        {
            y4[ m ] = y1[ m ] + h * k3[ m ];
        }

        BlasiusFun( y4, k4 );

        for ( int m = 0; m < nvec; ++ m )
        {
            y1[ m ] += coef * ( k1[ m ] + 2 * k2[ m ] + 2 * k3[ m ] + k4[ m ] );
        }
        eta += h;
        cout << eta << " " << y1[ 0 ] << " " << y1[ 1 ] << " " << y1[ 2 ] << "\n";
    }
    double alpha = pow( y1[ 1 ], - 3 / 2.0 );
    cout << "alpha = " << alpha << " nstep = " << n << "\n";
}


EndNameSpace