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

#include "CircleLineMesh.h"
#include "LineInfo.h"
#include "PointMachine.h"
#include "SegmentCtrl.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

CircleLineMesh::CircleLineMesh()
{
}

CircleLineMesh::~CircleLineMesh()
{
}

void CircleLineMesh::GenerateLineMesh()
{
    ;
}

void CircleLineMesh::CalcCurveGeometry()
{
    PointType * pt1 = point_Machine.GetPoint( this->curveInfo->p1 );
    PointType * pt2 = point_Machine.GetPoint( this->curveInfo->p2 );
    PointType * cp = point_Machine.GetPoint( this->center );

    Real x0 = pt1->x;
    Real y0 = pt1->y;
    Real z0 = pt1->z;

    Real x1 = pt2->x;
    Real y1 = pt2->y;
    Real z1 = pt2->z;

    Real xc = cp->x;
    Real yc = cp->y;
    Real zc = cp->z;

    Real dx0 = x0 - xc;
    Real dy0 = y0 - yc;

    Real dx1 = x1 - xc;
    Real dy1 = y1 - yc;

    Real r0 = DIST( dx0, dy0 );
    Real r1 = DIST( dx1, dy1 );

    this->alpha0 = acos( dx0 / r0 );
    if ( dy0 < 0 ) alpha0 = 2 * PI - alpha0;
    this->alpha1 = acos( dx1 / r1 );
    if ( dy1 < 0 ) alpha1 = 2 * PI - alpha1;

    if ( ( alpha0 - alpha1 ) > PI )
    {
        alpha1 += 2.0 * PI;
    }
    else if ( ( alpha1 - alpha0 ) > PI )
    {
        alpha1 -= 2.0 * PI;
    }

    this->radius = half * ( r0 + r1 );

    this->segmentCtrl->lenth = radius * ( alpha1 - alpha0 );
}

void CircleLineMesh::CalcCoor( Real s, Real & xt, Real & yt, Real & zt )
{
    PointType * pt1 = point_Machine.GetPoint( this->curveInfo->p1 );
    PointType * pt2 = point_Machine.GetPoint( this->curveInfo->p2 );
    PointType * cp  = point_Machine.GetPoint( this->center );

    Real angleSpan = ( this->alpha1 - this->alpha0 );
    Real ratio = s / this->segmentCtrl->lenth;

    Real cit = alpha0 + ratio * angleSpan;
    xt = cp->x + radius * cos( cit );
    yt = cp->y + radius * sin( cit );
    zt = half * ( pt1->z + pt2->z );
}

EndNameSpace