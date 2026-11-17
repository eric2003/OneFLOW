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

    solverInfo->stageSolutionNames.clear();
    for ( std::size_t i = 1; i < unsteady.flow.size(); ++ i )
    {
        solverInfo->stageSolutionNames.push_back( unsteady.flow[ i ] );
    }

    solverInfo->stageResidualNames.clear();
    for ( std::size_t i = 1; i < unsteady.residual.size(); ++ i )
    {
        solverInfo->stageResidualNames.push_back( unsteady.residual[ i ] );
    }
}

void RegisterSolverFieldInfo(
    SolverInfo * solverInfo,
    int nEqu,
    int nTEqu,
    const SolverLogisticsNames & logistics,
    const UnsteadyFieldNames & unsteady,
    const std::string & resFileName )
{
    solverInfo->nEqu  = nEqu;
    solverInfo->nTEqu = nTEqu;
    ApplyFieldNameContract( solverInfo, logistics, unsteady );
    solverInfo->resFileName = resFileName;
}

EndNameSpace
