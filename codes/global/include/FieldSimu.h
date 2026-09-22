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

BeginNameSpace( ONEFLOW )

class SimuContext;

// ------------------------------------------------------------
// Field-solve stages (fixed order; no numerical meaning change).
// ------------------------------------------------------------
void FieldSimuSetupGlobals();
void FieldSimuLoadGrid();
void FieldSimuPrepareWallDist();
void FieldSimuCreateSolvers();
void FieldSimuCreateSolvers( const SimuContext & ctx );
void FieldSimuInitFlowField();
void DumpFieldEnvironments();
void DumpCommunicationEnvironments();
void FieldSimuRun();

// Explicit pipeline: one place documents stage order and runs them.
// Free functions FieldSimuRunPipeline / FieldSimu remain as thin wrappers
// so existing call sites (SolveFieldTask, legacy) need not change.
struct FieldPipeline
{
    static constexpr int kStageCount = 6;

    // Human-readable stage tags (source form; not MessageMap names).
    static constexpr const char * StageName( int index )
    {
        switch ( index )
        {
        case 0: return "SetupGlobals";
        case 1: return "LoadGrid";
        case 2: return "PrepareWallDist";
        case 3: return "CreateSolvers";
        case 4: return "InitFlowField";
        case 5: return "Run";
        default: return "";
        }
    }

    // Production entry (no SimuContext injection).
    static void Run();

    // Context-aware entry: CreateSolvers uses expanded solver names when set.
    static void Run( SimuContext & ctx );
};

// Thin wrappers (same order as FieldPipeline::Run).
void FieldSimuRunPipeline();
void FieldSimuRunPipeline( SimuContext & ctx );

// Convenience: RunPipeline only (legacy / non-registry callers).
void FieldSimu();

void InitFlowSimuGlobal();
void InitializeSolver();

EndNameSpace
