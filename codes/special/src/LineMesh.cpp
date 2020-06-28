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

#include "LineMesh.h"
#include "LineInfo.h"
#include "SegmentCtrl.h"
#include "PointMachine.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

LineMesh::LineMesh()
{
}

LineMesh::~LineMesh()
{
}

void LineMesh::GenerateLineMesh()
{
    if ( ! this->IsValidState() ) return;
    this->GenerateCurveMesh();
}

void LineMesh::CalcCurveGeometry()
{
    PointType * pt1 = point_Machine.GetPoint( this->curveInfo->p1 );
    PointType * pt2 = point_Machine.GetPoint( this->curveInfo->p2 );

    Real x0 = pt1->x;
    Real y0 = pt1->y;
    Real z0 = pt1->z;

    Real x1 = pt2->x;
    Real y1 = pt2->y;
    Real z1 = pt2->z;

    Real dx = ( x1 - x0 );
    Real dy = ( y1 - y0 );
    Real dz = ( z1 - z0 );

    this->segmentCtrl->lenth = DIST( dx, dy, dz );
}

void LineMesh::CalcCoor( Real s, Real & xt, Real & yt, Real & zt )
{
    PointType * pt1 = 0;
    PointType * pt2 = 0;

    if ( this->segmentCtrl->c1 != 0 )
    {
        pt1 = point_Machine.GetPoint( this->curveInfo->p1 );
        pt2 = point_Machine.GetPoint( this->curveInfo->p2 );
    }
    else
    {
        pt1 = point_Machine.GetPoint( this->curveInfo->p2 );
        pt2 = point_Machine.GetPoint( this->curveInfo->p1 );
    }

    Real x0 = pt1->x;
    Real y0 = pt1->y;
    Real z0 = pt1->z;

    Real x1 = pt2->x;
    Real y1 = pt2->y;
    Real z1 = pt2->z;

    Real dx = ( x1 - x0 );
    Real dy = ( y1 - y0 );
    Real dz = ( z1 - z0 );

    Real ratio = s / this->segmentCtrl->lenth;

    xt = pt1->x + ratio * dx;
    yt = pt1->y + ratio * dy;
    zt = pt1->z + ratio * dz;
}

EndNameSpace