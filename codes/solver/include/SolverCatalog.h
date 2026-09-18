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
#include "SolverMap.h"

BeginNameSpace( ONEFLOW )

// Semantic facade over the process-default solver directory.
//
// Today ownership still lives in SolverMap (unique_ptr buckets + index maps).
// Call sites that mean "the simulation's solver catalog" should prefer this
// name so a later move to Session/SimuContext-owned storage is a type rename
// rather than a hunt for SolverMap:: statics.
//
// GetSolver returns a non-owning view (same as SolverMap::GetSolver).
struct SolverCatalog
{
    static void CreateDefault()
    {
        SolverMap::CreateSolvers();
    }

    static void CreateDefault( int gridType, const StringField * names )
    {
            SolverMap::CreateSolvers( gridType, names );
    }
    
    static void FreeDefault()
    {
        SolverMap::FreeSolverMap();
    }
    
    static Solver * GetSolver( int solverIndex, int gridType )
    {
        return SolverMap::GetSolver( solverIndex, gridType );
    }
};

EndNameSpace
