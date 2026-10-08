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

#include "GridGeneration.h"
#include "GridConversion.h"
#include "DomainInp.h"
#include "ClassicGrid.h"
#include "Partition.h"
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
    void GenerateClassic(
        const GridConfig & config,
        const std::string & caseDir )
    {
        GenerateClassicGrid( config );
        ConvertGrid( config, caseDir );
    }

    void ConvertOnly(
        const GridConfig & config,
        const std::string & caseDir )
    {
        ConvertGrid( config, caseDir );
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

    // The public functions only select a workflow. Concrete workflow steps do
    // not depend on a GridGeneration instance.
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
