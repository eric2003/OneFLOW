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

#include "INsSolver.h"
#include "INsCom.h"
#include "INsCtrl.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "UINsBcSolver.h"
#include "HeatFlux.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( INsSolver )

bool INsSolver::initFlag = false;

INsSolver::INsSolver()
{
}

INsSolver::~INsSolver()
{
}

void INsSolver::StaticInit()
{
    if ( INsSolver::initFlag ) return;
    INsSolver::initFlag = true;
    this->sTid = ONEFLOW::INC_NS_SOLVER;

    ins_ctrl.Init();
    inscom.Init();

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->sTid );
    solverInfo->nEqu  = inscom.nEqu;
    solverInfo->nTEqu = inscom.nTEqu;
    solverInfo->registerInterface = 0;
    solverInfo->residualName = "res";
    solverInfo->resFileName = GetDataValue< string >( "resFile" );
    solverInfo->gradString.push_back( "q"    );
    solverInfo->gradString.push_back( "dqdx" );
    solverInfo->gradString.push_back( "dqdy" );
    solverInfo->gradString.push_back( "dqdz" );

    solverInfo->implicitString.push_back( "q"  );
    solverInfo->implicitString.push_back( "dq" );
}

EndNameSpace