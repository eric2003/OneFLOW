/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\\*---------------------------------------------------------------------------*/

#include "SimuContext.h"

#include "Prj.h"
#include "SimuContextBootstrap.h"

BeginNameSpace( ONEFLOW )

SimuContext::SimuContext( const std::string & caseDir, bool debug )
    : caseDir_( Prj::ResolveCaseDir( caseDir ) )
{
    Prj::hx_debug = debug;
    Prj::run_from_ide = debug;
}

void SimuContext::SetupCaseEnvironment()
{
    // Readiness means the complete case configuration has been applied.
    // A failed setup is still cleaned up by the owning RunCase/rollback path.
    envReady_ = false;

    SimuContextDetail::InitializeCaseConfiguration( caseDir_ );

    // Publish readiness only after the full bootstrap sequence succeeds.
    envReady_ = true;
}

EndNameSpace
