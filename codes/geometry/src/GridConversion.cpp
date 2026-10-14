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

#include "GridConversion.h"
#include "CgnsFactory.h"
#include "GridMediator.h"
#include "Su2Grid.h"
#include "Plot3D.h"
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
}

namespace
{
bool IsSupportedConversion( GridFileType sourceType, GridFileType targetType )
{
    switch ( sourceType )
    {
        case GridFileType::Plot3D:
            return targetType == GridFileType::OneFLOW ||
                   targetType == GridFileType::CGNS;
        case GridFileType::SU2:
            return targetType == GridFileType::OneFLOW;
        case GridFileType::CGNS:
            return targetType == GridFileType::OneFLOW ||
                   targetType == GridFileType::CGNS;
        default:
            return false;
    }
}
} // namespace

void ConvertGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    if ( ! IsSupportedConversion( config.sourceType, config.targetType ) )
    {
        throw std::invalid_argument(
            std::string( "Unsupported grid conversion: " ) +
            std::string( ToString( config.sourceType ) ) + " -> " +
            std::string( ToString( config.targetType ) ) );
    }
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

EndNameSpace
