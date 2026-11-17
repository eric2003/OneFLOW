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

#pragma once

#include "NamespaceMacros.h"
#include <string>

BeginNameSpace( ONEFLOW )

namespace SimuContextDetail
{

// Execute the production case bootstrap against the legacy runtime boundary.
void InitializeCaseConfiguration( const std::string & caseDir );

}

EndNameSpace
