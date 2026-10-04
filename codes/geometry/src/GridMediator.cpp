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
#include "GridTypes.h"
#include "DataBase.h"
#include <utility>

BeginNameSpace( ONEFLOW )

void GridMediator::ReadGrid()
{
    // Dispatch by typed format (case-insensitive via ParseGridFileType).
    switch ( ParseGridFileType( this->gridType ) )
    {
        case GridFileType::Gridgen:
            this->ReadGridgen();
            break;
        case GridFileType::Plot3D:
            this->ReadPlot3D();
            break;
        default:
            // Unsupported or empty type: no-op (historical behavior).
            break;
    }
}

void GridMediator::AddDefaultName()
{
    const int numberOfZones = this->numberOfZones;

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( &GridAt( this->gridVector, iZone ) );

        grid->name = AddString( "Zone", iZone + 1 );

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup.get();
        const int nBcRegions = static_cast< int >( bcRegionGroup->regions.size() );
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

void ZgridMediator::add( std::unique_ptr< GridMediator > mediator )
{
    if ( ! mediator )
    {
        throw std::invalid_argument( "ZgridMediator cannot own a null GridMediator" );
    }

    this->mediators_.push_back( std::move( mediator ) );
}

GridMediator & ZgridMediator::at( int index )
{
    return *this->mediators_.at( static_cast< size_t >( index ) );
}

const GridMediator & ZgridMediator::at( int index ) const
{
    return *this->mediators_.at( static_cast< size_t >( index ) );
}

int ZgridMediator::size() const noexcept
{
    return static_cast< int >( this->mediators_.size() );
}

bool ZgridMediator::empty() const noexcept
{
    return this->mediators_.empty();
}

std::string ZgridMediator::targetFile() const
{
    return this->mediators_.at( 0 )->targetFile;
}

void ZgridMediator::CreateSimple( int nZone )
{
    auto gridMediator = std::make_unique< GridMediator >();
    gridMediator->numberOfZones = nZone;
    this->add( std::move( gridMediator ) );
}

void ZgridMediator::ReadGrid()
{
    this->ReadGrid( GridConfig::FromDataBase() );
}

void ZgridMediator::ReadGrid( const GridConfig & config )
{
    auto gridMediator = std::make_unique< GridMediator >();
    gridMediator->gridFile = config.sourceFile;
    gridMediator->bcFile   = config.bcFile;
    gridMediator->gridType = std::string( ToString( config.sourceType ) );

    // sourceCaseDir identifies where the input grid and its boundary file live.
    // An empty value keeps the historical current-project behavior.
    gridMediator->caseDir = config.sourceCaseDir;

    gridMediator->ReadGrid();
    this->add( std::move( gridMediator ) );
}

EndNameSpace
