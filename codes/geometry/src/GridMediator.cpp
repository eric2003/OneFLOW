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

#include "GridMediator.h"
#include "Plot3D.h"
#include "Su2Grid.h"
#include "StrGrid.h"
#include "StrUtil.h"
#include "BcRecord.h"
#include "GridPara.h"


using namespace std;

BeginNameSpace( ONEFLOW )


GridMediator::GridMediator()
{
    ;
}

GridMediator::~GridMediator()
{
    ;
}

void GridMediator::ReadGrid()
{
    if ( this->gridType == "gridgen" )
    {
        this->ReadGridgen();
    }
    else if ( this->gridType == "plot3d" )
    {
        this->ReadPlot3D();
    }
}

void GridMediator::AddDefaultName()
{
    int numberOfZones = this->numberOfZones;

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( this->gridVector[ iZone ] );

        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        grid->name = AddString( "Zone", iZone + 1 );

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        int nBcRegions = bcRegionGroup->regions->size();
        int icount = 0;
        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );
            bcRegion->regionName = AddString( "R", ir + 1 );

            int bcType = bcRegion->bcType;
            if ( bcType < 0 )
            {
                int sz = iZone + 1;
                int tz = bcRegion->t->zid + 1;
                bcRegion->regionName = AddString( "I", icount + 1 );
                icount ++;
            }
        }
    }
}

void GridMediator::ReadPlot3D()
{
    Plot3D::ReadPlot3D( this );
}

void GridMediator::ReadPlot3DCoor()
{
    Plot3D::ReadCoor( this );
}

void GridMediator::ReadGridgen()
{
}

ZgridMediator::ZgridMediator()
{
    this->flag = false;
}

ZgridMediator::~ZgridMediator()
{
    if ( this->flag )
    {
        for ( int i = 0; i < this->gm.size(); ++ i )
        {
            delete this->gm[ i ];
        }
    }
}

void ZgridMediator::AddGridMediator( GridMediator * gridMediator )
{
    this->gm.push_back( gridMediator );
}

GridMediator * ZgridMediator::GetGridMediator( int iGridMediator )
{
    return this->gm[ iGridMediator ];
}

int ZgridMediator::GetSize()
{
    return this->gm.size();
}

void ZgridMediator::CreateSimple( int nZone )
{
    GridMediator * gridMediator = new GridMediator();
    gridMediator->numberOfZones = nZone;
    this->AddGridMediator( gridMediator );
    this->flag = true;
}

void ZgridMediator::ReadGrid()
{
    GridMediator * gridMediator = new GridMediator();
    gridMediator->gridFile = grid_para.gridFile;
    gridMediator->bcFile = grid_para.bcFile;

    gridMediator->gridType = grid_para.filetype;
    gridMediator->ReadGrid();
    this->AddGridMediator( gridMediator );
    this->flag = true;
}

string ZgridMediator::GetTargetFile()
{
    int index = 0;
    GridMediator * gridMediator = this->gm[ index ];
    return gridMediator->targetFile;
}

void ZgridMediator::SetDeleteFlag( bool flag )
{
    this->flag = flag;
}

GridMediator * GlobalGrid::gridMediator = 0;

GlobalGrid::GlobalGrid()
{
    ;
}

GlobalGrid::~GlobalGrid()
{
    ;
}

void GlobalGrid::SetCurrentGridMediator( GridMediator * gridMediatorIn )
{
    GlobalGrid::gridMediator = gridMediatorIn;
}

Grid * GlobalGrid::GetGrid( int zoneId )
{
    return GlobalGrid::gridMediator->gridVector[ zoneId ];
}

EndNameSpace