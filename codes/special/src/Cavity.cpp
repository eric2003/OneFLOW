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

#include "Cavity.h"
#include "CurveLine.h"
#include "CurveMesh.h"
#include "BasicIO.h"
#include "Prj.h"
#include "DataBaseIO.h"
#include "Boundary.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Cavity::Cavity()
{
    ;
}

Cavity::~Cavity()
{
    ;
}

void Cavity::Run()
{
    int ni = 101;
    int nj = 51;
    int nNode = ni * nj;

    RealField2D x;
    RealField2D y;

    AllocateVector( x, ni, nj );
    AllocateVector( y, ni, nj );

    Real xl = 0.0;
    Real xr = 1.0;
    Real yl = 0.0;
    Real yr = 1.0;

    Real dx = ( xr - xl ) / ( ni - 1 );
    Real dy = ( yr - yl ) / ( nj - 1 );

    RealField xN;
    RealField yN;
    RealField zN;
    //Generate an evenly spaced, planar grid.
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            x[ i ][ j ] = xl + i * dx;
            y[ i ][ j ] = yl + j * dy;
            xN.push_back( x[ i ][ j ] );
            yN.push_back( y[ i ][ j ] );
            zN.push_back( 0 );
        }
    }

    fstream file;
    OpenPrjFile( file, "/grid/cavity2d.grd", ios_base::out|ios_base::binary );
    int nZone = 1;
    int nk = 1;
    HXWrite( & file, nZone );
    HXWrite( & file, ni );
    HXWrite( & file, nj );
    HXWrite( & file, nk );

    HXWrite( & file, xN );
    HXWrite( & file, yN );
    HXWrite( & file, zN );

    CloseFile( file );

    OpenPrjFile( file, "/grid/cavity2d.inp", ios_base::out );
    int solver = 1;
    string zName = "A";
    int nBc = 4;
    file << solver << endl;
    file << nZone << endl;
    file << ni << " " << nj << endl;
    file << zName << endl;
    file << nBc << endl;
    file << 1 << " " << ni  << " " << 1  << " " << 1  << " " << BC::SOLID_SURFACE << endl;
    file << 1  << " " << ni  << " " << nj  << " " << nj  << " " << BC::SOLID_SURFACE << endl;
    file << 1  << " " << 1  << " " << 1  << " " << nj  << " " << BC::SOLID_SURFACE << endl;
    file << ni  << " " << ni  << " " << 1  << " " << nj  << " " << BC::SOLID_SURFACE << endl;
    CloseFile( file );
}


EndNameSpace