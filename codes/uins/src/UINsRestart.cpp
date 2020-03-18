/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "UINsRestart.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "Iteration.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "INsCom.h"
#include "Ctrl.h"
#include "CellMesh.h"

BeginNameSpace( ONEFLOW )

UINsRestart::UINsRestart()
{
    ;
}

UINsRestart::~UINsRestart()
{
    ;
}

void UINsRestart::InitinsRestart( int sTid )
{
    Iteration::outerSteps = 0;
    ctrl.currTime = 0.0;

    UnsGrid * grid = Zone::GetUnsGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, "q" );

    int nEqu = inscom.inflow.size();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        SetField( ( * q )[ iEqu ], inscom.inflow[ iEqu ] );
    }
    this->InitUnsteady( sTid );

    RwInterface( sTid, GREAT_ZERO );

    if ( ctrl.inflowType == 3 )
    {
        RealField & xcc = grid->cellMesh->xcc;
        RealField & ycc = grid->cellMesh->ycc;
        RealField & zcc = grid->cellMesh->zcc;

        Real a = ctrl.initplane[ 0 ];
        Real b = ctrl.initplane[ 1 ];
        Real c = ctrl.initplane[ 2 ];
        Real d = ctrl.initplane[ 3 ];

        int nTCell = xcc.size();

        for ( int cId = 0; cId < nTCell; ++ cId )
        {
            Real x = xcc[ cId ];
            Real y = ycc[ cId ];
            Real z = zcc[ cId ];
            Real s = a * x + b * y + c * z + d;
            if ( s < 0 )
            {
                for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
                {
                    ( * q )[ iEqu ][ cId ] = ctrl.initflow1[ iEqu ];
                }
            }
            else
            {
                for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
                {
                    ( * q )[ iEqu ][ cId ] = ctrl.initflow2[ iEqu ];
                }
            }
        }
    }
    else if ( ctrl.inflowType == 4 )
    {
        RealField & xcc = grid->cellMesh->xcc;
        RealField & ycc = grid->cellMesh->ycc;
        RealField & zcc = grid->cellMesh->zcc;

        int nTCell = xcc.size();

        for ( int cId = 0; cId < nTCell; ++ cId )
        {
            Real x = xcc[ cId ];
            Real y = ycc[ cId ];
            Real z = zcc[ cId ];
        }
    }
}




EndNameSpace