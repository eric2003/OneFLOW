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

#include "CurveMesh.h"
#include "LineInfo.h"
#include "SegmentCtrl.h"
#include "PointMachine.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

CurveMesh::CurveMesh()
{
    this->state = 0;
}

CurveMesh::~CurveMesh()
{
    for ( int i = 0; i < segmentCtrl->nPoint; ++ i )
    {
        delete ptList[ i ];
    }
}

bool CurveMesh::IsValidState()
{
    if ( this->state == 1 ) return false;

    if ( segmentCtrl->distribution == 2 )
    {
    }
    return true;
}

PointType & CurveMesh::GetPoint( int id, int signFlag )
{
    int coef = ( 1 + signFlag ) / 2;
    int index1 = id;
    int index2 = ptList.size() - 1 - index1;
    int index = coef * index1 + ( 1 - coef ) * index2;
    return * ptList[ index ];
}

void CurveMesh::GenerateCurveMesh()
{
    ptList.resize( segmentCtrl->nPoint );
    for ( int i = 0; i < segmentCtrl->nPoint; ++ i )
    {
        ptList[ i ] = new PointType();
    }

    int st = 0;
    int ed = segmentCtrl->nPoint - 1;

    PointType * pt1 = point_Machine.GetPoint( this->curveInfo->p1 );
    PointType * pt2 = point_Machine.GetPoint( this->curveInfo->p2 );

    * ptList[ st ] = * pt1;
    * ptList[ ed ] = * pt2;

    this->CalcCurveGeometry();
    segmentCtrl->CalcFactor();

    for ( int iPoint = 1; iPoint < segmentCtrl->nPoint - 1; ++ iPoint )
    {
        Real factor = this->segmentCtrl->factorList[ iPoint ];
        int idx = this->segmentCtrl->pidxList[ iPoint ];

        Real s = this->segmentCtrl->lenth * factor;

        Real xt, yt, zt;
        this->CalcCoor( s, xt, yt, zt );

        ptList[ idx ]->x = xt;
        ptList[ idx ]->y = yt;
        ptList[ idx ]->z = zt;
    }

    this->state = 1;
}

void Alloc( RealField2D & field, int ni, int nj )
{
    field.resize( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        field[ i ].resize( nj );
    }
}

EndNameSpace