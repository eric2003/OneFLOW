/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
using namespace std;

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

class PointBasic
{
public:
    PointBasic();
    ~PointBasic();
public:
    typedef Point< Real > PointType;
    typedef set< PointType, PointCompare< Real > > PointSet;
public:
    PointSet pointSet;
    HXVector< PointType > pointList;
public:
    UInt GetNPoint() { return pointList.size(); }
    int AddPoint( Real xm, Real ym, Real zm );
    int DeletePoint( Real xm, Real ym, Real zm );
    int DeletePoint( PointBasic::PointType & point );
    int FindPoint( Real xm, Real ym, Real zm );
    void GetFaceCoorList( IntField & nodeId, RealField &xList, RealField &yList, RealField &zList );
protected:
    int AddPoint( PointBasic::PointType & point );
    int FindPointId( PointBasic::PointType & point );
    bool FindPoint( PointBasic::PointType & point, PointBasic::PointSet::iterator & iter );
private:
    void AddPointDirectly( PointBasic::PointType & point, int bcType = 0 );
};

class PointFactory : public PointBasic
{
public:
    PointFactory();
    ~PointFactory();
public:
    IntField c2g;
public:
    void InitC2g();
};

EndNameSpace