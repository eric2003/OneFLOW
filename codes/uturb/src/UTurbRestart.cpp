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

#include "UTurbRestart.h"
#include "SolverDef.h"
#include "Iteration.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "Ctrl.h"
#include "CellMesh.h"

BeginNameSpace( ONEFLOW )

UTurbRestart::UTurbRestart()
{
    ;
}

UTurbRestart::~UTurbRestart()
{
    ;
}

void UTurbRestart::InitRestart( int sTid )
{
    Iteration::outerSteps = 0;
    ctrl.currTime = 0.0;

    UnsGrid * grid = Zone::GetUnsGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, "turbq" );

    int nEqu = turbcom.inflow.size();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        SetField( ( * q )[ iEqu ], turbcom.inflow[ iEqu ] );
    }
    this->InitUnsteady( sTid );

    RwInterface( sTid, GREAT_ZERO );
}


EndNameSpace