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

#pragma once
#include "HXDefine.h"
#include "UnsteadyFieldNames.h"
#include <string>

BeginNameSpace( ONEFLOW )

// ---------------------------------------------------------------------------
// Pure field-name contracts (no Zone / SolverState / DataBase).
// Intent is encoded here so unit tests can lock production naming without a
// full simulation bootstrap. Keep SolverInfo::StaticInit and
// system/*/alloc/unsteady.txt in sync with Expected* helpers.
// ---------------------------------------------------------------------------

// Dual-time history order: Current, Previous, Old.
//   flow[0]      = solution          (q / turbq)
//   flow[1..2]   = stage solutions   (q1,q2 / turbq1,turbq2)
//   residual[0]  = residual          (res / turbres)
//   residual[1..2]= stage residuals  (res1,res2 / turbres1,turbres2)
UnsteadyFieldNames MakeUnsteadyFieldNames(
    const std::string & solution,
    const StringField & stageSolutions,
    const std::string & residual,
    const StringField & stageResiduals );

UnsteadyFieldNames ExpectedNsUnsteadyFieldNames();
UnsteadyFieldNames ExpectedINsUnsteadyFieldNames();
UnsteadyFieldNames ExpectedTurbUnsteadyFieldNames();

// ---------------------------------------------------------------------------
// Solver logistics name table (mirrors SolverInfo string slots).
// solutionGradNames: [0]=solution, [1..3]=cell-center grad dx,dy,dz
// implicitFieldNames: [0]=solution, [1]=update delta
// bcGradFieldNames: [0..2]=BC face grad dx,dy,dz
// ---------------------------------------------------------------------------
struct SolverLogisticsNames
{
    StringField solutionGradNames;
    StringField implicitFieldNames;
    StringField bcGradFieldNames;
    std::string residualName;
    std::string limiterName;
    std::string bcSolutionName;
};

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
    const std::string & bcSolution );

SolverLogisticsNames ExpectedNsSolverLogisticsNames();
SolverLogisticsNames ExpectedINsSolverLogisticsNames();
SolverLogisticsNames ExpectedTurbSolverLogisticsNames();

EndNameSpace
