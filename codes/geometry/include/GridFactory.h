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
#include "GridTypes.h"
#include <string>

BeginNameSpace( ONEFLOW )

// Offline grid generation / conversion / partition entry.
// Workflow execution is handled by the dispatch tables in GridFactory.cpp.
void GenerateGrid();

// Load the legacy database configuration and run it for a specific case.
void GenerateGrid( const std::string & caseDir );

// Run an explicit configuration.
void GenerateGrid( const GridConfig & config );

// Run an explicit configuration for a specific case.
void GenerateGrid( const GridConfig & config, const std::string & caseDir );

EndNameSpace
