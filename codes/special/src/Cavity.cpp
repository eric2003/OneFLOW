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
#include "BcRecord.h"
#include "Dimension.h"
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

    Real xl = 0.0;
    Real xr = 1.0;
    Real yl = 0.0;
    Real yr = 1.0;

    Real dx = ( xr - xl ) / ( ni - 1 );
    Real dy = ( yr - yl ) / ( nj - 1 );

    //Generate an evenly spaced, planar grid.
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
                Real xx = xl + ii * dx;
                Real yy = yl + jj * dy;
                Real zz = 0.0;

                xs( i, j, k ) = xx;
                ys( i, j, k ) = yy;
                zs( i, j, k ) = zz;
            }
        }
    }

    BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
    int nBcRegions = 4;
    grid->bcRegionGroup->Create( nBcRegions );

    BcRegion* bcRegion = 0;
    int ir = 0;
    bcRegion = new BcRegion( iZone, ir );
    bcRegion->s->SetRegion( 1, ni, 1, 1 );
    bcRegion->s->zid = iZone;
    bcRegion->regionName = "lower";
    bcRegion->bcType = BC::SOLID_SURFACE;
    bcRegionGroup->SetBcRegion( ir, bcRegion );
    ++ ir;

    bcRegion = new BcRegion( iZone, ir );
    bcRegion->s->SetRegion( 1, ni, nj, nj );
    bcRegion->s->zid = iZone;
    bcRegion->regionName = "upper";
    bcRegion->bcType = BC::SOLID_SURFACE;
    bcRegionGroup->SetBcRegion( ir, bcRegion );
    ++ ir;

    bcRegion = new BcRegion( iZone, ir );
    bcRegion->s->SetRegion( 1, 1, 1, nj );
    bcRegion->s->zid = iZone;
    bcRegion->regionName = "left";
    bcRegion->bcType = BC::SOLID_SURFACE;
    bcRegionGroup->SetBcRegion( ir, bcRegion );
    ++ ir;

    bcRegion = new BcRegion( iZone, ir );
    bcRegion->s->SetRegion( ni, ni, 1, nj );
    bcRegion->s->zid = iZone;
    bcRegion->regionName = "right";
    bcRegion->bcType = BC::SOLID_SURFACE;
    bcRegionGroup->SetBcRegion( ir, bcRegion );
    ++ ir;

    this->DumpPlot3DGrid( gridMediator );

    this->DumpCgnsGrid( gridMediator );

    delete gridMediator;
}

void Cavity::DumpPlot3DGrid( GridMediator * gridMediator )
{
    fstream file;
    OpenPrjFile( file, "/grid/cavity2d.grd", ios_base::out|ios_base::binary );
    int nZone = gridMediator->numberOfZones;
    HXWrite( & file, nZone );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * gridIn = gridMediator->gridVector[ iZone ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridIn );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;
        HXWrite( & file, ni );
        HXWrite( & file, nj );
        HXWrite( & file, nk );

        HXWrite( & file, grid->nodeMesh->xN );
        HXWrite( & file, grid->nodeMesh->yN );
        HXWrite( & file, grid->nodeMesh->zN );
    }

    CloseFile( file );

    OpenPrjFile( file, "/grid/cavity2d.inp", ios_base::out );
    int solver = 1;
    file << solver << "\n";
    file << nZone << "\n";
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * gridIn = gridMediator->gridVector[ iZone ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridIn );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        file << ni << " " << nj;
        if ( ONEFLOW::IsThreeD() )
        {
            file << nk;
        }
        file << "\n";
        file << grid->name << "\n";
        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        int nBcRegions = bcRegionGroup->regions->size();

        file << nBcRegions << "\n";

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );
            int bcType = bcRegion->bcType;
            BasicRegion * s = bcRegion->s;
            int imin = s->start[ 0 ];
            int imax = s->end[ 0 ];
            int jmin = s->start[ 1 ];
            int jmax = s->end[ 1 ];
            int kmin = s->start[ 2 ];
            int kmax = s->end[ 2 ];
            file << imin << " " << imax << " " << jmin << " " << jmax << " ";
            if ( ONEFLOW::IsThreeD() )
            {
                file << kmin << " " << kmax << " ";
            }
            file << bcType << "\n";
            if ( bcType < 0 )
            {
                BasicRegion * t = bcRegion->t;
                int imin = t->start[ 0 ];
                int imax = t->end[ 0 ];
                int jmin = t->start[ 1 ];
                int jmax = t->end[ 1 ];
                int kmin = t->start[ 2 ];
                int kmax = t->end[ 2 ];
                file << imin << " " << imax << " " << jmin << " " << jmax << " ";
                if ( ONEFLOW::IsThreeD() )
                {
                    file << kmin << " " << kmax << " ";
                }
                file << t->zid << "\n";
            }
        }

    }
    CloseFile( file );
}

void Cavity::DumpCgnsGrid( GridMediator * gridMediator )
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->DumpCgnsGrid( gridMediator );

    delete cgnsFactory;
}


EndNameSpace