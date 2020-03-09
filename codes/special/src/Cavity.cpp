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
#include "Plot3D.h"
#include "DataBase.h"
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
    gridMediator->gridFile = ONEFLOW::GetDataValue< string >( "sourceGridFileName" );
    gridMediator->bcFile   = ONEFLOW::GetDataValue< string >( "sourceGridBcName" );
    gridMediator->targetFile = ONEFLOW::GetDataValue< string >( "targetGridFileName" );

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

    BcRegion * bcRegion = 0;
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
    //bcRegion->bcType = BC::INTERFACE;
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
    Plot3D::DumpCoor( gridMediator );
    Plot3D::DumpBc( gridMediator );

}

void Cavity::DumpCgnsGrid( GridMediator * gridMediator )
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    ZgridMediator zgridMediator;
    zgridMediator.AddGridMediator( gridMediator );

    cgnsFactory->DumpCgnsGrid( & zgridMediator );

    delete cgnsFactory;
}


EndNameSpace