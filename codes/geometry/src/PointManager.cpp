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

#include "PointManager.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

PointManager::PointManager()
{
    // Use PointCompare's static tolerance, or wrap another layer here
}

PointManager::~PointManager()
{
}

void PointManager::SetTolerance(Real tol)
{
    PointCompare<Real>::ResetTolerance(tol);
}

PointManager::PointType PointManager::MakeKey(Real x, Real y, Real z) const
{
    // Deliberately set id to 0; we do not use Point's built-in id at all
    return PointType(x, y, z, 0);
}

int PointManager::AddPoint(Real x, Real y, Real z)
{
    PointType key = MakeKey(x, y, z);

    auto it = pointToId_.find(key);
    if (it != pointToId_.end())
    {
        return it->second;          // already exists
    }

    // new point
    int newId = static_cast<int>(points_.size());
    points_.push_back(key);         // The stored Point's id member is still 0; we don't care
    pointToId_[key] = newId;

    return newId;
}

int PointManager::FindPoint(Real x, Real y, Real z) const
{
    PointType key = MakeKey(x, y, z);
    auto it = pointToId_.find(key);
    if (it == pointToId_.end())
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return it->second;
}

void PointManager::GetPoint(int id, Real& x, Real& y, Real& z) const
{
    const PointType& pt = points_.at(id);   // at() performs bounds checking
    x = pt.x;
    y = pt.y;
    z = pt.z;
}

const PointManager::PointType& PointManager::GetPoint(int id) const
{
    return points_.at(id);
}

void PointManager::Clear()
{
    pointToId_.clear();
    points_.clear();
}

bool PointManager::DeletePoint(Real x, Real y, Real z)
{
    PointType key = MakeKey(x, y, z);
    auto it = pointToId_.find(key);
    if (it == pointToId_.end())
    {
        return false;
    }

    int id = it->second;
    pointToId_.erase(it);

    // Simple strategy: do not compact points_, just mark the slot as invalid (or leave a hole)
    // Choose the simplest approach here: set coordinates to an obviously invalid value, or keep as-is
    // (let the upper layer decide whether to reuse the id)
    // A more rigorous approach would be to maintain a free-list; keep it simple for now.
    return true;
}

bool PointManager::DeletePoint(int id)
{
    if (id < 0 || id >= static_cast<int>(points_.size()))
    {
        return false;
    }

    const PointType& pt = points_[id];
    PointType key = MakeKey(pt.x, pt.y, pt.z);

    auto it = pointToId_.find(key);
    if (it != pointToId_.end() && it->second == id)
    {
        pointToId_.erase(it);
        // Also do not compact the vector, leave a hole
        return true;
    }
    return false;
}

void PointManager::GetFaceCoorList( IntField & nodeIds, RealField &xList, RealField &yList, RealField &zList )
{
    xList.clear();
    yList.clear();
    zList.clear();
    xList.reserve(nodeIds.size());
    yList.reserve(nodeIds.size());
    zList.reserve(nodeIds.size());

    for (int ip : nodeIds)
    {
        const PointType& pt = points_[ip];
        xList.push_back(pt.x);
        yList.push_back(pt.y);
        zList.push_back(pt.z);
    }
}


MeshPointManager::MeshPointManager()
{
}

MeshPointManager::~MeshPointManager()
{
}

void MeshPointManager::InitLocalToGlobal()
{
    int nPoint = this->GetNPoint();
    this->localToGlobal.resize(nPoint);

    for (int iNode = 0; iNode < nPoint; ++iNode)
    {
        // Currently a simple identity mapping.
        // Later this can be customized when a grid only uses a subset of points.
        this->localToGlobal[iNode] = iNode;
    }
}

EndNameSpace
