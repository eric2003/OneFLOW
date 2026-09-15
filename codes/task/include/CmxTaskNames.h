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
#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

// CmxTask / MessageMap operation names.
// Must match action registration tables (script under system action/).
// Single source for FieldSimu, Multigrid, and contract tests.

inline constexpr const char* kInitFlowFieldTaskName   = "INIT_FLOWFIELD";
inline constexpr const char* kPostProcessTaskName     = "POST_PROCESS";
inline constexpr const char* kStoreRhsTaskName        = "STORE_RHS";
inline constexpr const char* kRestrictAllQTaskName    = "RESTRICT_ALL_Q";
inline constexpr const char* kLoadQTaskName           = "LOAD_Q";
inline constexpr const char* kLoadResidualsTaskName   = "LOAD_RESIDUALS";
inline constexpr const char* kUpdateResidualsTaskName = "UPDATE_RESIDUALS";
inline constexpr const char* kRestrictDefectTaskName  = "RESTRICT_DEFECT";
inline constexpr const char* kModifyCoarseGridTaskName  = "MODIFY_COARSEGRID";
inline constexpr const char* kModifyFineGridTaskName    = "MODIFY_FINEGRID";
inline constexpr const char* kRecoverCoarseGridTaskName = "RECOVER_COARSEGRID";
inline constexpr const char* kRecoverResidualsTaskName  = "RECOVER_RESIDUALS";
inline constexpr const char* kZeroResidualsTaskName     = "ZERO_RESIDUALS";

EndNameSpace