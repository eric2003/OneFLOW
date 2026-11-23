/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once
#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

// Optional resident macro-step driver for the validated device chain.
//
// The regular path dispatches one command-list task per stage phase
// (CALC_TIME_STEP, then UPDATE_RESIDUALS / CALC_LHS / UPDATE_FLOWFIELD /
// CALC_BOUNDARY per stage).  Every dispatch re-expands the solver function map
// into commands and resolves its parameters, which is host work that the
// 1D standalone port does not have: there the whole macro-step is one call.
//
// This driver reproduces the same sequence by calling the same entry points
// directly, so the numerics are unchanged, and lets us measure how much of the
// per-step time that dispatch layer costs.  Opt in with
// ONEFLOW_HIP_RESIDENT_MACROSTEP=1.
bool NsResidentMacroStepSupported();

// Advances exactly one outer step (all RK stages) through the device chain.
void NsResidentMacroStep();

EndNameSpace
