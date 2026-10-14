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
#include "SolverInfo.h"
#include <string>

BeginNameSpace( ONEFLOW )

// ---------------------------------------------------------------------------
// Generic field-name machinery only (no NS / INS / Turb knowledge).
//
// Each solver owns its Expected* tables next to StaticInit so solvers stay
// "happily ignorant" of each other. Cross-solver reads (e.g. turb -> mean-flow)
// use explicit solverType pins on FieldHome, not shared name tables.
// ---------------------------------------------------------------------------

// Dual-time history order: Current, Previous, Old.
UnsteadyFieldNames MakeUnsteadyFieldNames(
    const std::string & solution,
    const StringField & stageSolutions,
    const std::string & residual,
    const StringField & stageResiduals );

// Logistics name table (mirrors SolverInfo string slots).
// solutionGradNames: [0]=solution, [1..3]=cell-center grad
// implicitFieldNames: [0]=solution, [1]=update delta
// bcGradFieldNames: [0..2]=BC face grad
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

// Fill SolverInfo field-name slots. Does not modify resFileName.
void ApplyFieldNameContract(
    SolverInfo * solverInfo,
    const SolverLogisticsNames & logistics,
    const UnsteadyFieldNames & unsteady );

// Register field-name tables + equation counts onto SolverInfo.
// Pure w.r.t. Zone/DataBase: caller supplies nEqu and resFileName.
// Separates "metadata registration" from StaticInit resource loading (BC files, etc.).
void RegisterSolverFieldInfo(
    SolverInfo * solverInfo,
    int nEqu,
    int nTEqu,
    const SolverLogisticsNames & logistics,
    const UnsteadyFieldNames & unsteady,
    const std::string & resFileName );

EndNameSpace
