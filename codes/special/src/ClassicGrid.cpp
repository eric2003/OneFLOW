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
#include "GridCreate.h"
#include "DataBase.h"
#include "DataBaseIO.h"
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
    using GridGenerator = void ( * )();

    void RunCavity()
    {
        Cavity cavity;
        cavity.Run();
    }

    void RunRae2822()
    {
        Rae2822 rae2822;
        rae2822.Run();
    }

    void RunCylinder()
    {
        Cylinder cylinder;
        cylinder.Run( 3 );
    }

    void RunGridCreate()
    {
        GridCreate gridCreate;
        gridCreate.Run( 4 );
    }

    void RunCgnsTest()
    {
        CgnsTest cgnsTest;
        cgnsTest.Run();
    }

    struct GridGenerationEntry
    {
        int id;
        GridGenerator run;
    };

    constexpr GridGenerationEntry kGridGenerationEntries[] = {
        { 1, &RunCavity },
        { 2, &RunRae2822 },
        { 3, &RunCylinder },
        { 4, &RunGridCreate },
        { 5, &RunCgnsTest },
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

void ClassicGrid::Run() const
{
    const int generationId = GetDataValue< int >( "igene" );

    for ( const auto & entry : kGridGenerationEntries )
    {
        if ( entry.id == generationId )
        {
            entry.run();
            return;
        }
    }
}

EndNameSpace
