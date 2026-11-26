/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software; you can redistribute it and/or modify it
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

#include "StlMeshIO.h"
#include "StlSurfaceTopology.h"

#include <array>
#include <vector>

namespace ONEFLOW {

struct StlCutPlane
{
    std::array<double, 3> point{};
    std::array<double, 3> normal{1.0, 0.0, 0.0};
    // Retain the positive signed-distance half-space when true.
    bool keepPositive{true};
};

struct StlSurfaceCutOptions
{
    StlCutPlane inlet;
    StlCutPlane outlet{ {}, {1.0, 0.0, 0.0}, false };
    double tolerance{1.0e-8};
    bool capBoundaries{true};
};

struct StlSurfaceCutResult
{
    StlMesh wall;
    StlMesh inletCap;
    StlMesh outletCap;
    StlMesh closedSurface;
    std::vector<std::array<double, 3>> inletLoop;
    std::vector<std::array<double, 3>> outletLoop;
    StlSurfaceTopologyReport topology;
};

// Clip an arbitrary closed, connected, outward-oriented triangle shell by two
// planes. The source is treated only as triangles; no analytic primitive data
// is inspected. The current algorithm intentionally rejects branched or
// multiple cut loops rather than silently producing an invalid cap.
StlSurfaceCutResult CutStlSurfaceBetweenPlanes(
    const StlMesh& source,
    const StlSurfaceCutOptions& options);

} // namespace ONEFLOW
