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

#include "CurveMachine.h"
#include "PointMachine.h"
#include "CurveLine.h"
#include "LineInfo.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

CurveMachine curve_Machine;

CurveMachine::CurveMachine()
{
}

CurveMachine::~CurveMachine()
{
    for ( int i = 0; i < curveList.size(); ++ i )
    {
        delete curveList[ i ];
    }
}

void CurveMachine::AddLine( int id1, int id2 )
{
    CurveLine * curveLine = new CurveLine();
    curveLine->lineType = LINE;
    PointType * p1 = point_Machine.GetPoint( id1 );
    PointType * p2 = point_Machine.GetPoint( id2 );
    curveLine->start_p = * p1;
    curveLine->end_p   = * p2;
    curveList.push_back( curveLine );
}

void CurveMachine::AddCircle( int id1, int id2, int id3 )
{
    CurveLine * curveLine = new CurveLine();
    curveLine->lineType = CIRCLE;
    PointType * p1 = point_Machine.GetPoint( id1 );
    PointType * p2 = point_Machine.GetPoint( id2 );
    PointType * p3 = point_Machine.GetPoint( id3 );
    curveLine->start_p = * p1;
    curveLine->end_p = * p2;
    curveLine->center_p = * p3;
    curveList.push_back( curveLine );
}

void CurveMachine::AddParabolic( int id1, int id2 )
{
    CurveLine * curveLine = new CurveLine();
    curveLine->lineType = PARABOLIC;
    PointType * p1 = point_Machine.GetPoint( id1 );
    PointType * p2 = point_Machine.GetPoint( id2 );
    curveLine->start_p = * p1;
    curveLine->end_p = * p2;
    curveList.push_back( curveLine );
}

CurveLine *  CurveMachine::GetCurve( int curveId )
{
    return this->curveList[ curveId ];
}

EndNameSpace