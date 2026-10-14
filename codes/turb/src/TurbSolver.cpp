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

#include "TurbSolver.h"
#include "TurbBcSolver.h"
#include "BcData.h"
#include "SolverInfo.h"
#include "FieldNameContract.h"
#include "SolverDef.h"
#include "TurbCtrl.h"
#include "TurbCom.h"
#include "DataBase.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( TurbSolver )

bool TurbSolver::initFlag = false;

TurbSolver::TurbSolver()
{
}

TurbSolver::~TurbSolver()
{
}

void TurbSolver::Reset()
{
    TurbSolver::initFlag = false;
}

void TurbSolver::StaticInit()
{
    if ( TurbSolver::initFlag ) return;
    TurbSolver::initFlag = true;

    std::string fileName = "grid/turb_bc.txt";
    turb_bc_data.Init( fileName );

    turb_ctrl.Init();
    turbcom.Init();

    this->solverType = ONEFLOW::TURB_SOLVER;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->solverType );

    solverInfo->nEqu  = turbcom.nEqu;
    solverInfo->nTEqu = turbcom.nTEqu;

    // Field names from pure contract (FieldNameContract); single source of truth.
    const SolverLogisticsNames logistics = ExpectedTurbSolverLogisticsNames();
    const UnsteadyFieldNames unsteady = ExpectedTurbUnsteadyFieldNames();
    solverInfo->residualName = logistics.residualName;
    solverInfo->solutionGradNames = logistics.solutionGradNames;
    solverInfo->implicitFieldNames = logistics.implicitFieldNames;
    solverInfo->bcGradFieldNames = logistics.bcGradFieldNames;
    solverInfo->limiterName = logistics.limiterName;
    solverInfo->bcSolutionName = logistics.bcSolutionName;
    solverInfo->resFileName = GetDataValue< std::string >( "turbresFile" );
    solverInfo->stageSolutionNames.clear();
    solverInfo->stageSolutionNames.push_back( unsteady.flow[ 1 ] );
    solverInfo->stageSolutionNames.push_back( unsteady.flow[ 2 ] );
    solverInfo->stageResidualNames.clear();
    solverInfo->stageResidualNames.push_back( unsteady.residual[ 1 ] );
    solverInfo->stageResidualNames.push_back( unsteady.residual[ 2 ] );
}

EndNameSpace
