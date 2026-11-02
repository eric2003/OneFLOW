/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "GridGeneration.h"
#include "GridWorkflow.h"
#include "ClassicGridWorkflow.h"
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
    struct PipelineEntry
    {
        GridObjective objective;
        void ( * run )(
            const GridConfig &,
            const std::string & );
    };

    // The dispatcher selects an explicit workflow implementation.
    constexpr PipelineEntry kPipelines[] = {
        { GridObjective::GenerateClassic, &ExecuteClassicGridWorkflow },
        { GridObjective::ConvertOnly,     &ExecuteGridConversion },
        { GridObjective::GenerateInp,     &ExecuteDomainInpWorkflow },
        { GridObjective::Partition,       &ExecutePartitionWorkflow },
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

void GenerateGrid()
{
    GenerateGrid( GridConfig::FromDataBase() );
}

void GenerateGrid( const std::string & caseDir )
{
    GenerateGrid( GridConfig::FromDataBase(), caseDir );
}

void GenerateGrid( const GridConfig & config )
{
    DispatchPipeline( config, "" );
}

void GenerateGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    DispatchPipeline( config, caseDir );
}

EndNameSpace
