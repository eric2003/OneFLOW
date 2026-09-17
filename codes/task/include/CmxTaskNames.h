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

// --- already present (keep) ---
// kInitFlowFieldTaskName, kPostProcessTaskName, kStoreRhsTaskName,
// kRestrictAllQTaskName, kLoadQTaskName, kLoadResidualsTaskName,
// kUpdateResidualsTaskName, kRestrictDefectTaskName,
// kModifyCoarseGridTaskName, kModifyFineGridTaskName,
// kRecoverCoarseGridTaskName, kRecoverResidualsTaskName,
// kZeroResidualsTaskName,

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

// TimeIntegral / LUSGS
inline constexpr const char* kCalcTimeStepTaskName         = "CALC_TIME_STEP";
inline constexpr const char* kCalcLhsTaskName              = "CALC_LHS";
inline constexpr const char* kUpdateFlowFieldTaskName      = "UPDATE_FLOWFIELD";
inline constexpr const char* kCalcBoundaryTaskName         = "CALC_BOUNDARY";
inline constexpr const char* kZeroDqFieldTaskName          = "ZERO_DQ_FIELD";
inline constexpr const char* kInitLusgsTaskName            = "INIT_LUSGS";
inline constexpr const char* kLusgsLowerSweepTaskName      = "LUSGS_LOWER_SWEEP";
inline constexpr const char* kExchangeInterfaceDqTaskName  = "EXCHANGE_INTERFACE_DQ";
inline constexpr const char* kLusgsUpperSweepTaskName      = "LUSGS_UPPER_SWEEP";
inline constexpr const char* kUpdateFlowFieldLusgsTaskName = "UPDATE_FLOWFIELD_LUSGS";
inline constexpr const char* kSolTurbTaskName              = "SOL_TURB";
inline constexpr const char* kSolHeatTaskName              = "SOL_HEAT";

// SolverState
inline constexpr const char* kCalcUnsteadyCriterionTaskName = "CALC_UNSTEADY_CRITERION";

// MultiBlock / wall distance
inline constexpr const char* kCalcMetricsTaskName       = "CALC_METRICS";
inline constexpr const char* kFillWallStructTaskName    = "FILL_WALL_STRUCT";
inline constexpr const char* kCalcWallDistTaskName      = "CALC_WALL_DIST";
inline constexpr const char* kWriteWallDistTaskName     = "WRITE_WALL_DIST";
inline constexpr const char* kReadWallDistTaskName      = "READ_WALL_DIST";
inline constexpr const char* kAllocateWallDistTaskName  = "ALLOCATE_WALL_DIST";

// Restart / init-flowfield command list (AddCmdToList path)
inline constexpr const char* kInitFirstTaskName      = "INIT_FIRST";
inline constexpr const char* kInitRestartTaskName    = "INIT_RESTART";
inline constexpr const char* kReadRestartTaskName    = "READ_RESTART";
inline constexpr const char* kInitInsRestartTaskName = "INIT_INSRESTART";
inline constexpr const char* kReadInsRestartTaskName = "READ_INSRESTART";
inline constexpr const char* kInitFinalTaskName      = "INIT_FINAL";

// Interface exchange
inline constexpr const char* kUploadInterfaceDataTaskName   = "UPLOAD_INTERFACE_DATA";
inline constexpr const char* kUpdateInterfaceDataTaskName   = "UPDATE_INTERFACE_DATA";
inline constexpr const char* kDownloadInterfaceDataTaskName = "DOWNLOAD_INTERFACE_DATA";

// Dump / visualization / unsteady (Ns / INs / Turb SolverImp)
inline constexpr const char* kDumpResidualTaskName       = "DUMP_RESIDUAL";
inline constexpr const char* kDumpAerodynamicTaskName    = "DUMP_AERODYNAMIC";
inline constexpr const char* kDumpPressureCoeffTaskName  = "DUMP_PRESSURE_COEFF";
inline constexpr const char* kDumpHeatfluxCoeffTaskName  = "DUMP_HEATFLUX_COEFF";
inline constexpr const char* kDumpRestartTaskName        = "DUMP_RESTART";
inline constexpr const char* kDumpLaminarPlateTaskName   = "DUMP_LAMINAR_PLATE";
inline constexpr const char* kDumpTurbPlateTaskName      = "DUMP_TURB_PLATE";
inline constexpr const char* kVisualizationTaskName      = "VISUALIZATION";
inline constexpr const char* kUpdateUnsteadyFlowTaskName = "UPDATE_UNSTEADY_FLOW";

EndNameSpace