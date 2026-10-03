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

#include "ClassicGrid.h"
#include "GridTypes.h"
#include "GridCreate.h"
#include "Boundary.h"
#include "HXMath.h"
#include "Cavity.h"
#include "Rae2822.h"
#include "Cylinder.h"
#include "CgnsTest.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

namespace
{
    using GridGenerator = void ( * )( const GridConfig & );

    void RunCavity( const GridConfig & )
    {
        Cavity cavity;
        cavity.Run();
    }

    void RunRae2822( const GridConfig & )
    {
        Rae2822 rae2822;
        rae2822.Run();
    }

    void RunCylinder( const GridConfig & )
    {
        Cylinder cylinder;
        cylinder.Run();
    }

    void RunGridCreate( const GridConfig & config )
    {
        GridCreate gridCreate;
        gridCreate.Run( config );
    }

    void RunCgnsTest( const GridConfig & )
    {
        CgnsTest cgnsTest;
        cgnsTest.Run();
    }

    struct GridGenerationEntry
    {
        GridGenerationType type;
        GridGenerator run;
    };

    constexpr GridGenerationEntry kGridGenerationEntries[] = {
        { GridGenerationType::Cavity,     &RunCavity },
        { GridGenerationType::Rae2822,    &RunRae2822 },
        { GridGenerationType::Cylinder,   &RunCylinder },
        { GridGenerationType::GridCreate, &RunGridCreate },
        { GridGenerationType::CgnsTest,   &RunCgnsTest },
    };
}

ClassicGrid::ClassicGrid()
{
    ;
}

ClassicGrid::~ClassicGrid()
{
    ;
}

void ClassicGrid::Run( const GridConfig & config ) const
{
    if ( ! config.generationType )
    {
        return;
    }

    for ( const auto & entry : kGridGenerationEntries )
    {
        if ( entry.type == *config.generationType )
        {
            entry.run( config );
            return;
        }
    }
}

EndNameSpace
