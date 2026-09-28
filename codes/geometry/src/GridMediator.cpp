/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "StringUtils.h"
#include "BcRecord.h"
#include "GridPara.h"
#include <utility>

BeginNameSpace( ONEFLOW )

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
    const int numberOfZones = this->numberOfZones;

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( this->gridVector[ iZone ] );

        grid->name = AddString( "Zone", iZone + 1 );

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        const int nBcRegions = static_cast< int >( bcRegionGroup->regions->size() );
        int icount = 0;
        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );
            bcRegion->regionName = AddString( "R", ir + 1 );

            const int bcType = bcRegion->bcType;
            if ( bcType < 0 )
            {
                bcRegion->regionName = AddString( "I", icount + 1 );
                ++ icount;
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

void ZgridMediator::AddGridMediator( GridMediator * gridMediator )
{
    this->gm.emplace_back( gridMediator );
}

void ZgridMediator::AddGridMediator( std::unique_ptr< GridMediator > gridMediator )
{
    this->gm.push_back( std::move( gridMediator ) );
}

GridMediator * ZgridMediator::GetGridMediator( int iGridMediator ) const
{
    return this->gm[ static_cast< size_t >( iGridMediator ) ].get();
}

int ZgridMediator::GetSize() const
{
    return static_cast< int >( this->gm.size() );
}

void ZgridMediator::CreateSimple( int nZone )
{
    auto gridMediator = std::make_unique< GridMediator >();
    gridMediator->numberOfZones = nZone;
    this->AddGridMediator( std::move( gridMediator ) );
}

void ZgridMediator::ReadGrid()
{
    auto gridMediator = std::make_unique< GridMediator >();
    gridMediator->gridFile = grid_para.gridFile;
    gridMediator->bcFile   = grid_para.bcFile;
    gridMediator->gridType = grid_para.filetype;
    gridMediator->ReadGrid();
    this->AddGridMediator( std::move( gridMediator ) );
}

std::string ZgridMediator::GetTargetFile() const
{
    const int index = 0;
    return this->gm[ static_cast< size_t >( index ) ]->targetFile;
}

GridMediator * GlobalGrid::gridMediator = nullptr;

void GlobalGrid::SetCurrentGridMediator( GridMediator * gridMediatorIn )
{
    GlobalGrid::gridMediator = gridMediatorIn;
}

GridMediator * GlobalGrid::GetCurrentGridMediator() noexcept
{
    return GlobalGrid::gridMediator;
}

Grid * GlobalGrid::GetGrid( int zoneId )
{
    GridMediator * gm = GlobalGrid::gridMediator;
    if ( ! gm )
    {
        throw std::logic_error(
            "GlobalGrid::GetGrid: no current GridMediator "
            "(call SetCurrentGridMediator or ScopedCurrentGridMediator first)" );
    }
    return gm->gridVector[ zoneId ];
}

EndNameSpace
