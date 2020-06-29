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

#include "ULhs.h"
#include "Zone.h"
#include "Ctrl.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "UNsCom.h"
#include "UCom.h"
#include "FieldWrap.h"
#include "SolverInfo.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ULhs::ULhs()
{
    ;
}

ULhs::~ULhs()
{
    ;
}

void ULhs::CalcLHS( int sTid )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );
    string & residualName = solverInfo->residualName;

    MRField * dq = FieldHome::GetUnsField( grid, residualName );

    for ( int iEqu = 0; iEqu < solverInfo->nTEqu; ++ iEqu )
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            //Real dt = ctrl.pdt;
            Real dt = ( * unsf.timestep )[ 0 ][ cId ];
            ( * dq )[ iEqu ][ cId ] *= ( dt / ( * ug.cvol )[ cId ] ) * ctrl.lhscoef;
        }
    }
}

EndNameSpace