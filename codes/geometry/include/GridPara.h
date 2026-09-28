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

// Mesh parameters.
// New code should prefer GridConfig::FromDataBase(). This class remains as a
// thin compatibility facade that mirrors the historical global grid_para.
class GridPara
{
public:
    GridPara() = default;
    ~GridPara() = default;

public:
    // Historical string fields (kept for call sites that still read them).
    std::string topo;
    std::string filetype;        // source format token
    std::string target_filetype; // target format token
    std::string format;
    std::string gridFile;
    std::string bcFile;
    std::string targetFile;

    // Historical integer objective (gridObj). Prefer objective() below.
    int gridObj{ 0 };

    int multiBlock{ 0 };
    int axis_dir{ 0 };
    Real gridScale{ 1.0 };
    RealField gridTrans;

public:
    // Load from DataBase and refresh both typed and legacy fields.
    void Init();

    // Typed accessors (C++20-friendly API).
    [[nodiscard]] GridObjective objective() const noexcept
    {
        return ParseGridObjective( gridObj ).value_or( GridObjective::ConvertOnly );
    }

    [[nodiscard]] GridFileType sourceType() const noexcept
    {
        return ParseGridFileType( filetype );
    }

    [[nodiscard]] GridFileType targetType() const noexcept
    {
        return ParseGridFileType( target_filetype );
    }

    // Build a GridConfig snapshot from the current (already Init'd) state.
    [[nodiscard]] GridConfig ToConfig() const;
};

extern GridPara grid_para;

int GetGridTopoType();

EndNameSpace
