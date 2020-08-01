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

#include "NsSolver.h"
#include "NsCom.h"
#include "NsCtrl.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "UNsBcSolver.h"
#include "HeatFlux.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( NsSolver )

bool NsSolver::initFlag = false;

NsSolver::NsSolver()
{
}

NsSolver::~NsSolver()
{
}

void NsSolver::StaticInit()
{
    if ( NsSolver::initFlag ) return;
    NsSolver::initFlag = true;
    this->sTid = ONEFLOW::NS_SOLVER;

    ns_ctrl.Init();
    nscom.Init();

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->sTid );
    solverInfo->nEqu  = nscom.nEqu;
    solverInfo->nTEqu = nscom.nTEqu;
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