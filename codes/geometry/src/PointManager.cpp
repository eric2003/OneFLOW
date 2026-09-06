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
    // 使用 PointCompare 的静态 tolerance，或在这里再包一层
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
    // 故意把 id 设为 0，我们完全不使用 Point 自带的 id
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
    points_.push_back(key);         // 这里存的 Point 的 id 成员仍是 0，我们不关心
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
    const PointType& pt = points_.at(id);   // at() 带边界检查
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

    // 简单策略：不压缩 points_，只是把该位置标记为无效（或直接留下空洞）
    // 这里选择最简单的做法：把坐标设成一个明显无效值，或保持原样（由上层决定是否复用 id）
    // 更严谨的做法是维护一个 free-list，这里先保持简单。
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
        // 同样不压缩 vector，留下空洞
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
        const PointType& pt = points_[ip];   // 注意用你现在的成员名
        xList.push_back(pt.x);
        yList.push_back(pt.y);
        zList.push_back(pt.z);
    }
}


PointFactory::PointFactory()
{
}

PointFactory::~PointFactory()
{
}

//void PointFactory::InitC2g()
//{
//    int nPoint = this->GetNPoint();   // use public interface instead of pointList.size()
//    this->c2g.resize(nPoint);
//
//    for (int iNode = 0; iNode < nPoint; ++iNode)
//    {
//        this->c2g[iNode] = iNode;
//    }
//}

void PointFactory::InitLocalToGlobal()
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
