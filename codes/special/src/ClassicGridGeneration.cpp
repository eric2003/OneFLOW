/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#include "ClassicGridGeneration.h"
#include "GridCreate.h"
#include "Cavity.h"
#include "Rae2822.h"
#include "Cylinder.h"
#include "CgnsTest.h"
#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace
{
    using GridGenerator = void ( * )( const GridConfig & );

    void RunCavity( const GridConfig & )
    {
        GenerateCavityGrid();
    }

    void RunRae2822( const GridConfig & )
    {
        GenerateRae2822Grid();
    }

    void RunCylinder( const GridConfig & )
    {
        Cylinder cylinder;
        cylinder.Run();
    }

    void RunGridCreate( const GridConfig & config )
    {
        GenerateLayoutGrid( config );
    }

    void RunCgnsTest( const GridConfig & )
    {
        CgnsTest cgnsTest;
        cgnsTest.Run();
    }

    enum class ClassicGeneratorId : int
    {
        Cavity     = 1,
        Rae2822    = 2,
        Cylinder   = 3,
        GridCreate = 4,
        CgnsTest   = 5
    };

    struct GridGenerationEntry
    {
        ClassicGeneratorId id;
        GridGenerator run;
    };

    constexpr GridGenerationEntry kGridGenerationEntries[] = {
        { ClassicGeneratorId::Cavity,     &RunCavity },
        { ClassicGeneratorId::Rae2822,    &RunRae2822 },
        { ClassicGeneratorId::Cylinder,   &RunCylinder },
        { ClassicGeneratorId::GridCreate, &RunGridCreate },
        { ClassicGeneratorId::CgnsTest,   &RunCgnsTest },
    };
}

void GenerateClassicGrid( const GridConfig & config )
{
    if ( ! config.generationId )
    {
        return;
    }

    const auto generatorId =
        static_cast< ClassicGeneratorId >( *config.generationId );

    for ( const auto & entry : kGridGenerationEntries )
    {
        if ( entry.id == generatorId )
        {
            entry.run( config );
            return;
        }
    }

    throw std::invalid_argument(
        "Unknown classic grid generation id: " +
        std::to_string( *config.generationId ) );
}

EndNameSpace
