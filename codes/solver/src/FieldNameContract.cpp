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

#include "FieldNameContract.h"

BeginNameSpace( ONEFLOW )

UnsteadyFieldNames MakeUnsteadyFieldNames(
    const std::string & solution,
    const StringField & stageSolutions,
    const std::string & residual,
    const StringField & stageResiduals )
{
    UnsteadyFieldNames names;
    names.flow.push_back( solution );
    for ( std::size_t i = 0; i < stageSolutions.size(); ++ i )
    {
        names.flow.push_back( stageSolutions[ i ] );
    }
    names.residual.push_back( residual );
    for ( std::size_t i = 0; i < stageResiduals.size(); ++ i )
    {
        names.residual.push_back( stageResiduals[ i ] );
    }
    return names;
}

UnsteadyFieldNames ExpectedNsUnsteadyFieldNames()
{
    StringField stageQ;
    stageQ.push_back( "q1" );
    stageQ.push_back( "q2" );
    StringField stageR;
    stageR.push_back( "res1" );
    stageR.push_back( "res2" );
    return MakeUnsteadyFieldNames( "q", stageQ, "res", stageR );
}

UnsteadyFieldNames ExpectedINsUnsteadyFieldNames()
{
    // INS shares NS dual-time names in system/ins/alloc/unsteady.txt.
    return ExpectedNsUnsteadyFieldNames();
}

UnsteadyFieldNames ExpectedTurbUnsteadyFieldNames()
{
    StringField stageQ;
    stageQ.push_back( "turbq1" );
    stageQ.push_back( "turbq2" );
    StringField stageR;
    stageR.push_back( "turbres1" );
    stageR.push_back( "turbres2" );
    return MakeUnsteadyFieldNames( "turbq", stageQ, "turbres", stageR );
}

SolverLogisticsNames MakeSolverLogisticsNames(
    const std::string & solution,
    const std::string & dqdx,
    const std::string & dqdy,
    const std::string & dqdz,
    const std::string & updateDelta,
    const std::string & residual,
    const std::string & bcdqdx,
    const std::string & bcdqdy,
    const std::string & bcdqdz,
    const std::string & limiter,
    const std::string & bcSolution )
{
    SolverLogisticsNames names;
    names.solutionGradNames.push_back( solution );
    names.solutionGradNames.push_back( dqdx );
    names.solutionGradNames.push_back( dqdy );
    names.solutionGradNames.push_back( dqdz );
    names.implicitFieldNames.push_back( solution );
    names.implicitFieldNames.push_back( updateDelta );
    names.bcGradFieldNames.push_back( bcdqdx );
    names.bcGradFieldNames.push_back( bcdqdy );
    names.bcGradFieldNames.push_back( bcdqdz );
    names.residualName = residual;
    names.limiterName = limiter;
    names.bcSolutionName = bcSolution;
    return names;
}

SolverLogisticsNames ExpectedNsSolverLogisticsNames()
{
    return MakeSolverLogisticsNames(
        "q", "dqdx", "dqdy", "dqdz", "dq", "res",
        "bcdqdx", "bcdqdy", "bcdqdz",
        "limiter", "bc_q" );
}

SolverLogisticsNames ExpectedINsSolverLogisticsNames()
{
    // INS logistics names match NS in INsSolver::StaticInit.
    return ExpectedNsSolverLogisticsNames();
}

SolverLogisticsNames ExpectedTurbSolverLogisticsNames()
{
    return MakeSolverLogisticsNames(
        "turbq", "turbdqdx", "turbdqdy", "turbdqdz", "turbdq", "turbres",
        "turbbcdqdx", "turbbcdqdy", "turbbcdqdz",
        "turblimiter", "" );
}

void ApplyFieldNameContract(
    SolverInfo * solverInfo,
    const SolverLogisticsNames & logistics,
    const UnsteadyFieldNames & unsteady )
{
    solverInfo->residualName = logistics.residualName;
    solverInfo->solutionGradNames = logistics.solutionGradNames;
    solverInfo->implicitFieldNames = logistics.implicitFieldNames;
    solverInfo->bcGradFieldNames = logistics.bcGradFieldNames;
    solverInfo->limiterName = logistics.limiterName;
    solverInfo->bcSolutionName = logistics.bcSolutionName;

    // unsteady.flow[0] / residual[0] are current; stages are Previous/Old.
    solverInfo->stageSolutionNames.clear();
    if ( unsteady.flow.size() >= 3 )
    {
        solverInfo->stageSolutionNames.push_back( unsteady.flow[ 1 ] );
        solverInfo->stageSolutionNames.push_back( unsteady.flow[ 2 ] );
    }
    solverInfo->stageResidualNames.clear();
    if ( unsteady.residual.size() >= 3 )
    {
        solverInfo->stageResidualNames.push_back( unsteady.residual[ 1 ] );
        solverInfo->stageResidualNames.push_back( unsteady.residual[ 2 ] );
    }
}

EndNameSpace
