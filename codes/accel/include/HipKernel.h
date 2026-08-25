/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

namespace ONEFLOW
{

// Executes a small double-precision device calculation.  It is opt-in and is
// used to validate the complete HIP/DCU compilation and runtime path before a
// production CFD kernel is enabled.
void RunHipBackendSelfTest();

}
