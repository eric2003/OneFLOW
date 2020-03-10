/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "DataBase.h"
#include "DataBaseIO.h"
#include "Boundary.h"
#include "HXMath.h"
#include "Sod.h"
#include "Cavity.h"
#include "Rae2822.h"
#include "Cylinder.h"
#include "CgnsTest.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ClassicGrid::ClassicGrid()
{
    ;
}

ClassicGrid::~ClassicGrid()
{
    ;
}

void ClassicGrid::Run()
{
    int igene = GetDataValue< int >( "igene" );
    if ( igene == 0 )
    {
        Sod * sod = new Sod();
        sod->Run();
        delete sod;
    }
    else if ( igene == 1 )
    {
        Cavity * cavity = new Cavity();
        cavity->Run();
        delete cavity;
    }
    else if ( igene == 2 )
    {
        Rae2822 * rae2822 = new Rae2822();
        rae2822->Run();
        delete rae2822;
    }
    else if ( igene == 3 )
    {
        Cylinder * cylinder = new Cylinder();
        cylinder->Run( igene );
        delete cylinder;
    }
    else if ( igene == 4 )
    {
        Cylinder * cylinder = new Cylinder();
        cylinder->Run( igene );
        delete cylinder;
    }
    else if ( igene == 5 )
    {
        CgnsTest * cgnsTest = new CgnsTest();
        cgnsTest->Run();
        delete cgnsTest;
    }
}

EndNameSpace