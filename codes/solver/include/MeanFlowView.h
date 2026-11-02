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
#include "HXArray.h"
#include "SolverDef.h"

BeginNameSpace( ONEFLOW )

class Grid;

// ---------------------------------------------------------------------------
// Explicit cross-solver read of mean-flow (NS) fields.
//
// Intent: turbulence (and similar) may need density/velocity gradients from
// the NS solution without owning those fields. Binding always uses NS_SOLVER
// -- never ambient SolverState -- so turb stays loosely coupled.
//
// All pointers are non-owning views (FieldHome / database owns storage).
// ---------------------------------------------------------------------------
struct MeanFlowView
{
    MRField * q    = nullptr;
    MRField * dqdx = nullptr;
    MRField * dqdy = nullptr;
    MRField * dqdz = nullptr;

    // Contract: mean-flow logistics belong to NS_SOLVER.
    static constexpr int SolverType() { return NS_SOLVER; }

    // Bind non-owning views for the given grid.
    void Bind( Grid * grid );
};

EndNameSpace
