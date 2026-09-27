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

BeginNameSpace( ONEFLOW )

// Offline grid generation / conversion / partition entry.
// Dispatch is table-driven (see GridFactory.cpp); no magic switch on int.
class GridFactory
{
public:
    GridFactory() = default;
    ~GridFactory() = default;

    // Load config from DataBase and run the selected pipeline.
    void Run();

    // Run with an explicit config (preferred for tests and callers that
    // already hold a GridConfig).
    void Run( const GridConfig & config );

public:
    // Pipeline steps (also used as registry targets).
    void DataBaseGrid();
    void ConvertGrid( const GridConfig & config );
    void GeneInp();
    void PartGrid();

    // Format-specific convert helpers.
    void Plot3DProcess( const GridConfig & config );
    void SU2Process();
    void CGNSProcess();
};

// Public entry used by the rest of the code base.
void GenerateGrid();

EndNameSpace
