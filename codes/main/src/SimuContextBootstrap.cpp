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
----------------------------------------------------------------------------*/

#include "SimuContextBootstrap.h"

#include "SimuContextIO.h"
#include "Prj.h"
#include "ParaFile.h"
#include "LogFile.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

namespace SimuContextDetail
{

void InitializeCaseConfiguration( const std::string & caseDir )
{
    // Bind all legacy case-relative state before reading case-owned inputs.
    Prj::SetPrjBaseDir( caseDir );
    logFile.SetCaseDir( caseDir );
    ONEFLOW::ReadControlInfo( caseDir );

    // Validate the complete settings snapshot before changing effective roots.
    const CaseIOSettings settings = ReadCaseIOSettings(
        RequireGlobalDataBase().RequireDataPara() );
    ApplyCaseIOSettings( settings );
}

}

EndNameSpace
