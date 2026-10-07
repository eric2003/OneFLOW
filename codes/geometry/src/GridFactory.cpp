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
#include "Partition.h"
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
    using Converter = void ( * )(
        const GridConfig &,
        const std::string & );

    struct ConverterEntry
    {
        GridFileType sourceType;
        Converter convert;
    };

    void ConvertPlot3DToOneFLOW(
        const GridConfig & config,
        const std::string & /*caseDir*/ )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.CommonToOneFlowGrid( config );
    }

    void ConvertPlot3DToCGNS(
        const GridConfig & config,
        const std::string & caseDir )
    {
        CgnsFactory cgnsFactory;
        ZgridMediator zgridMediator;
        Plot3D::Plot3DToCgns( &zgridMediator, config, caseDir );
        cgnsFactory.DumpCgnsGrid( zgridMediator );
    }

    void ConvertPlot3D(
        const GridConfig & config,
        const std::string & caseDir )
    {
        switch ( config.targetType )
        {
            case GridFileType::OneFLOW:
                ConvertPlot3DToOneFLOW( config, caseDir );
                return;
            case GridFileType::CGNS:
                ConvertPlot3DToCGNS( config, caseDir );
                return;
            default:
                throw std::invalid_argument(
                    std::string( "Unsupported Plot3D target type: " ) +
                    std::string( ToString( config.targetType ) ) );
        }
    }

    void ConvertSU2(
        const GridConfig & config,
        const std::string & caseDir )
    {
        Su2Grid su2Grid;
        su2Grid.Su2ToOneFlowGrid( config, caseDir );
    }

    void ConvertCGNS(
        const GridConfig & config,
        const std::string & caseDir )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.GenerateGrid( config, caseDir );
    }

    constexpr ConverterEntry kConverters[] = {
        { GridFileType::Plot3D, &ConvertPlot3D },
        { GridFileType::SU2,    &ConvertSU2 },
        { GridFileType::CGNS,   &ConvertCGNS },
    };

    void DispatchConverter(
        const GridConfig & config,
        const std::string & caseDir )
    {
        for ( const auto & entry : kConverters )
        {
            if ( entry.sourceType == config.sourceType )
            {
                entry.convert( config, caseDir );
                return;
            }
        }

        throw std::invalid_argument(
            std::string( "Unsupported source grid type: " ) +
            std::string( ToString( config.sourceType ) ) );
    }

    void GenerateClassic(
        const GridConfig & config,
        const std::string & caseDir )
    {
        GenerateClassicGrid( config );
        DispatchConverter( config, caseDir );
    }

    void ConvertOnly(
        const GridConfig & config,
        const std::string & caseDir )
    {
        DispatchConverter( config, caseDir );
    }

    void GenerateInp(
        const GridConfig & /*config*/,
        const std::string & /*caseDir*/ )
    {
        DomainInp domainInp;
        domainInp.Run();
    }

    void PartitionGrid(
        const GridConfig & /*config*/,
        const std::string & /*caseDir*/ )
    {
        Partition part;
        part.Run();
    }

    struct PipelineEntry
    {
        GridObjective objective;
        void ( * run )(
            const GridConfig &,
            const std::string & );
    };

    // The factory only selects a workflow. Concrete workflow steps do not
    // depend on a GridFactory instance.
    constexpr PipelineEntry kPipelines[] = {
        { GridObjective::GenerateClassic, &GenerateClassic },
        { GridObjective::ConvertOnly,     &ConvertOnly },
        { GridObjective::GenerateInp,     &GenerateInp },
        { GridObjective::Partition,      &PartitionGrid },
    };

    void DispatchPipeline(
        const GridConfig & config,
        const std::string & caseDir )
    {
        for ( const auto & entry : kPipelines )
        {
            if ( entry.objective == config.objective )
            {
                entry.run( config, caseDir );
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
    GridFactory gf;
    gf.Run();
}

void GenerateGrid( const std::string & caseDir )
{
    GridFactory gf;
    gf.Run( GridConfig::FromDataBase(), caseDir );
}

void GridFactory::Run()
{
    Run( GridConfig::FromDataBase() );
}

void GridFactory::Run( const GridConfig & config )
{
    DispatchPipeline( config, "" );
}

void GridFactory::Run(
    const GridConfig & config,
    const std::string & caseDir )
{
    DispatchPipeline( config, caseDir );
}

EndNameSpace
