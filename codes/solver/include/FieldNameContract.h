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

// Pure field-name contracts (no Zone / SolverState / DataBase).
// Intent: dual-time history order is always Current, Previous, Old.
//
//   flow[0]      = solution          (q / turbq)
//   flow[1]      = stage solution 0  (q1 / turbq1)
//   flow[2]      = stage solution 1  (q2 / turbq2)
//   residual[0]  = residual          (res / turbres)
//   residual[1]  = stage residual 0  (res1 / turbres1)
//   residual[2]  = stage residual 1  (res2 / turbres2)
//
// Must stay aligned with system/*/alloc/unsteady.txt and SolverInfo StaticInit.

UnsteadyFieldNames MakeUnsteadyFieldNames(
    const std::string & solution,
    const StringField & stageSolutions,
    const std::string & residual,
    const StringField & stageResiduals );

// Documented production contracts (same strings as alloc/unsteady.txt).
UnsteadyFieldNames ExpectedNsUnsteadyFieldNames();
UnsteadyFieldNames ExpectedINsUnsteadyFieldNames();
UnsteadyFieldNames ExpectedTurbUnsteadyFieldNames();

EndNameSpace
