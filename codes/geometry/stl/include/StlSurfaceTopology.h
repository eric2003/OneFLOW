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

#include <cstddef>

namespace ONEFLOW {

struct StlSurfaceTopologyReport
{
    std::size_t weldedVertexCount{0};
    std::size_t edgeCount{0};
    std::size_t boundaryEdgeCount{0};
    std::size_t nonManifoldEdgeCount{0};
    std::size_t connectedComponentCount{0};
    double signedVolume{0.0};

    bool IsClosedTwoManifold() const noexcept
    {
        return boundaryEdgeCount == 0 && nonManifoldEdgeCount == 0;
    }

    bool HasConsistentOutwardOrientation(double tolerance = 0.0) const noexcept
    {
        return IsClosedTwoManifold() && signedVolume > tolerance;
    }
};

// Analyze a triangle shell using tolerance-based vertex welding. The input mesh
// is never modified. Disconnected components, open edges, non-manifold edges,
// and signed volume are reported independently.
StlSurfaceTopologyReport AnalyzeStlSurface(
    const StlMesh& mesh,
    double weldTolerance = 1.0e-9);

} // namespace ONEFLOW
