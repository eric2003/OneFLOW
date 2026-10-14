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

#include "NsSolver.h"
#include "NsCom.h"
#include "NsCtrl.h"
#include "SolverInfo.h"
#include "FieldNameContract.h"
#include "SolverDef.h"
#include "UNsBcSolver.h"
#include "HeatFlux.h"
#include "DataBase.h"
#include "BcData.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( NsSolver )

bool NsSolver::initFlag = false;

void NsSolver::Reset()
{
    NsSolver::initFlag = false;
}

void NsSolver::StaticInit()
{
    if ( NsSolver::initFlag ) return;
    NsSolver::initFlag = true;

    ns_ctrl.Init();
    nscom.Init();
    std::string fileName = "grid/ns_bc.txt";
    ns_bc_data.Init( fileName );

    this->solverType = ONEFLOW::NS_SOLVER;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->solverType );
    solverInfo->nEqu  = nscom.nEqu;
    solverInfo->nTEqu = nscom.nTEqu;
    // Field names: pure contract (unit-tested) -> SolverInfo.
    ApplyFieldNameContract(
        solverInfo,
        ExpectedNsSolverLogisticsNames(),
        ExpectedNsUnsteadyFieldNames() );
    solverInfo->resFileName = GetDataValue< std::string >( "resFile" );
}

EndNameSpace
