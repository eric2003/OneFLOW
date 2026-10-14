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

#include "INsSolver.h"
#include "INsBcSolver.h"
#include "BcData.h"
#include "INsCom.h"
#include "INsCtrl.h"
#include "SolverInfo.h"
#include "FieldNameContract.h"
#include "SolverDef.h"
#include "UINsBcSolver.h"
#include "HeatFlux.h"
#include "DataBase.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( INsSolver )

bool INsSolver::initFlag = false;

void INsSolver::Reset()
{
    INsSolver::initFlag = false;
}

void INsSolver::StaticInit()
{
    if ( INsSolver::initFlag ) return;
    INsSolver::initFlag = true;
    this->solverType = ONEFLOW::INC_NS_SOLVER;

    ins_ctrl.Init();
    nscom.Init();
    std::string fileName = "grid/ins_bc.txt";
    ins_bc_data.Init( fileName );

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->solverType );
    solverInfo->nEqu  = nscom.nEqu;
    solverInfo->nTEqu = nscom.nTEqu;

    // Field names: pure contract (unit-tested) -> SolverInfo.
    ApplyFieldNameContract(
        solverInfo,
        ExpectedINsSolverLogisticsNames(),
        ExpectedINsUnsteadyFieldNames() );
    solverInfo->resFileName = GetDataValue< std::string >( "resFile" );
}

EndNameSpace
