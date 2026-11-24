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

#include "StlSurfaceTopology.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace ONEFLOW {
namespace {

using Point = std::array<double, 3>;
using Cell = std::array<long long, 3>;
using Edge = std::pair<std::size_t, std::size_t>;

struct EdgeUse
{
    std::size_t count{0};
    int orientation{0};
};

double DistanceSquared(const Point& a, const Point& b)
{
    const double x = a[0] - b[0];
    const double y = a[1] - b[1];
    const double z = a[2] - b[2];
    return x*x + y*y + z*z;
}

Cell MakeCell(const Point& point, double tolerance)
{
    Cell cell{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const double scaled = std::floor(point[axis] / tolerance);
        if (!std::isfinite(scaled) ||
            scaled < static_cast<double>(std::numeric_limits<long long>::min() + 2) ||
            scaled > static_cast<double>(std::numeric_limits<long long>::max() - 2))
            throw std::invalid_argument("STL coordinate is outside the vertex-welding range");
        cell[axis] = static_cast<long long>(scaled);
    }
    return cell;
}

std::vector<std::size_t> Weld(
    const std::vector<Point>& points, double tolerance, std::size_t& uniqueCount)
{
    std::map<Cell, std::vector<std::size_t>> buckets;
    std::vector<Point> representatives;
    std::vector<std::size_t> ids;
    ids.reserve(points.size());
    const double toleranceSquared = tolerance * tolerance;

    for (const Point& point : points)
    {
        for (double value : point)
            if (!std::isfinite(value))
                throw std::invalid_argument("STL contains a non-finite vertex");

        const Cell cell = MakeCell(point, tolerance);
        std::size_t match = representatives.size();
        for (long long dx = -1; dx <= 1; ++dx)
        for (long long dy = -1; dy <= 1; ++dy)
        for (long long dz = -1; dz <= 1; ++dz)
        {
            const Cell neighbor{cell[0] + dx, cell[1] + dy, cell[2] + dz};
            const auto found = buckets.find(neighbor);
            if (found == buckets.end()) continue;
            for (std::size_t candidate : found->second)
            {
                if (DistanceSquared(point, representatives[candidate]) < toleranceSquared)
                    match = std::min(match, candidate);
            }
        }
        if (match == representatives.size())
        {
            match = representatives.size();
            representatives.push_back(point);
            buckets[cell].push_back(match);
        }
        ids.push_back(match);
    }
    uniqueCount = representatives.size();
    return ids;
}

} // namespace

StlSurfaceTopologyReport AnalyzeStlSurface(
    const StlMesh& mesh, double weldTolerance)
{
    if (!std::isfinite(weldTolerance) || weldTolerance <= 0.0)
        throw std::invalid_argument("STL weld tolerance must be positive and finite");
    if (mesh.triangles.empty())
        throw std::invalid_argument("Cannot analyze an STL mesh without triangles");

    std::vector<Point> vertices;
    vertices.reserve(mesh.triangles.size() * 3);
    double signedVolume = 0.0;
    for (const StlTriangle& triangle : mesh.triangles)
    {
        for (const Point& point : triangle.vertices) vertices.push_back(point);
        const Point& a = triangle.vertices[0];
        const Point& b = triangle.vertices[1];
        const Point& c = triangle.vertices[2];
        const Point ab{b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        const Point ac{c[0]-a[0], c[1]-a[1], c[2]-a[2]};
        const Point cross{
            ab[1]*ac[2]-ab[2]*ac[1],
            ab[2]*ac[0]-ab[0]*ac[2],
            ab[0]*ac[1]-ab[1]*ac[0]};
        const double twiceArea = std::sqrt(
            cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
        if (!std::isfinite(twiceArea) || twiceArea <= weldTolerance*weldTolerance)
            throw std::invalid_argument("STL contains a zero-area or near-degenerate triangle");
        signedVolume += (
            a[0] * (b[1]*c[2] - b[2]*c[1]) +
            a[1] * (b[2]*c[0] - b[0]*c[2]) +
            a[2] * (b[0]*c[1] - b[1]*c[0])) / 6.0;
    }

    StlSurfaceTopologyReport report;
    std::vector<std::size_t> ids = Weld(vertices, weldTolerance, report.weldedVertexCount);
    report.signedVolume = signedVolume;

    std::map<Edge, EdgeUse> edges;
    std::vector<std::vector<std::size_t>> adjacency(report.weldedVertexCount);
    for (std::size_t triangle = 0; triangle < mesh.triangles.size(); ++triangle)
    {
        const std::size_t base = triangle * 3;
        const std::array<std::size_t, 3> tri{ids[base], ids[base+1], ids[base+2]};
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
            throw std::invalid_argument("STL contains a degenerate triangle after vertex welding");
        for (std::size_t side = 0; side < 3; ++side)
        {
            const std::size_t a = tri[side];
            const std::size_t b = tri[(side + 1) % 3];
            const Edge key = std::minmax(a, b);
            EdgeUse& use = edges[key];
            ++use.count;
            use.orientation += (a < b) ? 1 : -1;
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }
    }

    report.edgeCount = edges.size();
    for (const auto& [edge, use] : edges)
    {
        (void)edge;
        if (use.count == 1) ++report.boundaryEdgeCount;
        else if (use.count != 2 || use.orientation != 0) ++report.nonManifoldEdgeCount;
    }

    std::vector<bool> visited(report.weldedVertexCount, false);
    for (std::size_t start = 0; start < visited.size(); ++start)
    {
        if (visited[start]) continue;
        ++report.connectedComponentCount;
        std::vector<std::size_t> pending{start};
        visited[start] = true;
        while (!pending.empty())
        {
            const std::size_t current = pending.back();
            pending.pop_back();
            for (std::size_t next : adjacency[current])
            {
                if (!visited[next])
                {
                    visited[next] = true;
                    pending.push_back(next);
                }
            }
        }
    }
    return report;
}

} // namespace ONEFLOW
