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

#include "GridFactory.h"
#include "CgnsFactory.h"
#include "GridMediator.h"
#include "DomainInp.h"
#include "Su2Grid.h"
#include "ClassicGrid.h"
#include "Plot3D.h"
#include "HXMath.h"
#include "Partition.h"
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
    // Free-function pipelines keep the registry simple (essence of the
    // suggested map + std::function design) while avoiding type-erasure
    // and heap cost: a plain function pointer table is enough.
    void PipelineGenerateClassic( GridFactory & self, const GridConfig & config )
    {
        self.DataBaseGrid();
        self.ConvertGrid( config );
    }

    void PipelineConvertOnly( GridFactory & self, const GridConfig & config )
    {
        self.ConvertGrid( config );
    }

    void PipelineGenerateInp( GridFactory & self, const GridConfig & /*config*/ )
    {
        self.GeneInp();
    }

    void PipelinePartition( GridFactory & self, const GridConfig & /*config*/ )
    {
        self.PartGrid();
    }

    struct PipelineEntry
    {
        GridObjective objective;
        void ( * run )( GridFactory &, const GridConfig & );
    };

    // Fixed-size, data-driven table. Easy to extend; no switch on magic int.
    constexpr PipelineEntry kPipelines[] = {
        { GridObjective::GenerateClassic, &PipelineGenerateClassic },
        { GridObjective::ConvertOnly,     &PipelineConvertOnly     },
        { GridObjective::GenerateInp,    &PipelineGenerateInp      },
        { GridObjective::Partition,      &PipelinePartition        },
    };

    void DispatchPipeline( GridFactory & self, const GridConfig & config )
    {
        for ( const auto & entry : kPipelines )
        {
            if ( entry.objective == config.objective )
            {
                entry.run( self, config );
                return;
            }
        }

        throw std::invalid_argument(
            std::string( "Unknown GridObjective / gridObj: " ) +
            std::string( ToString( config.objective ) ) );
    }
}

// Generates the grid based on the global configuration.
void GenerateGrid()
{
    // Stack allocation: no polymorphic need, automatic cleanup.
    GridFactory gf;
    gf.Run();
}

void GenerateGrid( const std::string & caseDir )
{
    // The case directory is explicit at the multi-case boundary.
    GridFactory gf;
    gf.Run( GridConfig::FromDataBase(), caseDir );
}

void GridFactory::Run()
{
    // Load the typed configuration directly from the database.
    Run( GridConfig::FromDataBase() );
}

void GridFactory::Run( const GridConfig & config )
{
    DispatchPipeline( *this, config );
}

void GridFactory::Run(
    const GridConfig & config,
    const std::string & caseDir )
{
    // Bind the case only for this grid operation; no global project state
    // is changed here.
    caseDir_ = caseDir;
    DispatchPipeline( *this, config );
}

void GridFactory::GeneInp()
{
    DomainInp domainInp;
    domainInp.Run();
}

void GridFactory::PartGrid()
{
    Partition part;
    part.Run();
}

void GridFactory::DataBaseGrid()
{
    ClassicGrid classicGrid;
    classicGrid.Run();
}

void GridFactory::ConvertGrid( const GridConfig & config )
{
    switch ( config.sourceType )
    {
        case GridFileType::Plot3D:
            this->Plot3DProcess( config );
            break;
        case GridFileType::SU2:
            this->SU2Process();
            break;
        case GridFileType::CGNS:
            this->CGNSProcess();
            break;
        default:
            throw std::invalid_argument(
                std::string( "Unsupported source grid type: " ) +
                std::string( ToString( config.sourceType ) ) );
    }
}


void GridFactory::Plot3DProcess( const GridConfig & config )
{
    if ( config.targetType == GridFileType::OneFLOW )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.SetCaseDir( caseDir_ );
        cgnsFactory.CommonToOneFlowGrid();
    }
    else if ( config.targetType == GridFileType::CGNS )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.SetCaseDir( caseDir_ );
        ZgridMediator zgridMediator;
        // Owned GridMediator instances are cleaned up automatically.
        Plot3D::Plot3DToCgns( &zgridMediator, caseDir_ );
        cgnsFactory.DumpCgnsGrid( zgridMediator );
    }
    else
    {
        throw std::invalid_argument(
            std::string( "Unsupported Plot3D target type: " ) +
            std::string( ToString( config.targetType ) ) );
    }
}

void GridFactory::SU2Process()
{
    Su2Grid su2Grid;
    su2Grid.SetCaseDir( caseDir_ );
    su2Grid.Su2ToOneFlowGrid();
}

void GridFactory::CGNSProcess()
{
    CgnsFactory cgnsFactory;
    cgnsFactory.SetCaseDir( caseDir_ );
    cgnsFactory.GenerateGrid();
}

EndNameSpace
