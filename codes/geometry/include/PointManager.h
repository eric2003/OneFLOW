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
#include "Point.h"
#include <vector>
#include <map>

BeginNameSpace( ONEFLOW )

template < typename T >
class PointCompare
{
public:
    typedef Point< T > point_type;
    static T tolerance;
    static void ResetTolerance( const T & toleranceIn );
public:
    bool operator()( const point_type & lhs, const point_type & rhs ) const
    {
        return lhs.Compare( rhs, PointCompare< T >::tolerance );
    }
};

template < typename T >
T PointCompare< T >::tolerance = static_cast< T >( 1.0e-10 );

template < typename T >
void PointCompare< T >::ResetTolerance( const T & toleranceIn )
{
    PointCompare< T >::tolerance = toleranceIn;
}

// --------------------------------------------------------------------------
// PointManager
// Manages unique 3D points with tolerance-based comparison.
// Important: We intentionally ignore Point::id. 
// Business IDs are managed externally by this class.
// --------------------------------------------------------------------------
class PointManager
{
public:
    using PointType = Point<Real>;          // 推荐用 PointType（更通用）
    // using PointKey  = Point<Real>;       // 如果想强调它是 map 的 key，也可以

    using PointCompareType = PointCompare<Real>;
    using PointMap         = std::map<PointType, int, PointCompareType>;  // geometry -> id
    using PointList        = std::vector<PointType>;                      // id -> geometry

public:
    PointManager();
    ~PointManager();

    // Reset tolerance used by PointCompare
    void SetTolerance(Real tol);

    // Add a point. Returns the unique business id.
    // If a point already exists within tolerance, returns the existing id.
    int AddPoint(Real x, Real y, Real z);

    // Find point. Returns id or -1 if not found.
    int FindPoint(Real x, Real y, Real z) const;

    // Get coordinates by business id
    void GetPoint(int id, Real& x, Real& y, Real& z) const;
    const PointType& GetPoint(int id) const;

    // Number of unique points
    int GetNPoint() const { return static_cast<int>(points_.size()); }

    // Optional: clear all
    void Clear();

    // Optional: delete (simple version, leaves id holes)
    // Returns true if the point was found and removed.
    bool DeletePoint(Real x, Real y, Real z);
    bool DeletePoint(int id);

    void GetFaceCoorList( IntField & nodeIds, RealField &xList, RealField &yList, RealField &zList );

protected:
    PointMap   pointToId_;   // geometry -> business id
    PointList  points_;      // business id -> geometry (Point::id is ignored)

    // Helper: create a temporary Point for query (id is set to 0, ignored)
    PointType MakeKey(Real x, Real y, Real z) const;
};

class PointFactory : public PointManager
{
public:
    PointFactory();
    ~PointFactory();

public:
    IntField localToGlobal;   // local node id -> global node id in PointManager

    void InitLocalToGlobal();
};
EndNameSpace
