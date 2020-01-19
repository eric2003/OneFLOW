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
#include "FileUtil.h"
#include "Prj.h"
#include "DataBaseIO.h"
#include "Boundary.h"
#include "HXMath.h"
#include "CgnsFactory.h"
#include "CgnsZone.h"
#include "GridMediator.h"
#include "BgGrid.h"
#include "StrGrid.h"
#include "StrUtil.h"
#include "NodeMesh.h"
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
    int nk = 1;

    int nZone = 1;
    GridMediator * gridMediator = new GridMediator();
    gridMediator->numberOfZones = nZone;
    gridMediator->gridVector.resize( nZone );

    Grid * gridstr = ONEFLOW::CreateStrGrid();
    StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
    int iZone = 0;
    gridMediator->gridVector[ iZone ] = grid;
    grid->name = AddString( "Zone", iZone );
    grid->id = iZone;
    grid->ni = ni;
    grid->nj = nj;
    grid->nk = nk;
    grid->SetBasicDimension();
    grid->nodeMesh->CreateNodes( grid->nNode );
    grid->SetLayout();

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

    Field3D & xs = * grid->strx;
    Field3D & ys = * grid->stry;
    Field3D & zs = * grid->strz;

    for ( int k = 1; k <= nk; ++ k )
    {
        for ( int j = 1; j <= nj; ++ j )
        {
            for ( int i = 1; i <= ni; ++ i )
            {
                int ii = i - 1;
                int jj = j - 1;
                Real xx = x[ ii ][ jj ];
                Real yy = y[ ii ][ jj ];
                Real zz = 0.0;

                xs( i, j, k ) = xx;
                ys( i, j, k ) = yy;
                zs( i, j, k ) = zz;
            }
        }
    }

    fstream file;
    OpenPrjFile( file, "/grid/cavity2d.grd", ios_base::out|ios_base::binary );
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

    this->DumpCgnsGrid( gridMediator );

    delete gridMediator;
}

void Cavity::DumpCgnsGrid( GridMediator * gridMediator )
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->DumpCgnsGrid( gridMediator );

    delete cgnsFactory;
}


EndNameSpace