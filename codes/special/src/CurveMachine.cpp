/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#include "CurveMachine.h"
#include "PointMachine.h"
#include "CurveLine.h"
#include "LineInfo.h"
#include "HXMath.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

CurveMachine curve_Machine;

CurveMachine::CurveMachine()
{
}

CurveMachine::~CurveMachine() = default;

void CurveMachine::AddLine( int id1, int id2 )
{
    auto curveLine = std::make_unique< CurveLine >();
    curveLine->lineType = LINE;
    const PointType & p1 = point_Machine.GetPoint( id1 );
    const PointType & p2 = point_Machine.GetPoint( id2 );
    curveLine->start_p = p1;
    curveLine->end_p   = p2;
    curveList.push_back( std::move( curveLine ) );
}

void CurveMachine::AddCircle( int id1, int id2, int id3 )
{
    auto curveLine = std::make_unique< CurveLine >();
    curveLine->lineType = CIRCLE;
    PointType * p1 = point_Machine.GetPoint( id1 );
    PointType * p2 = point_Machine.GetPoint( id2 );
    const PointType & p3 = point_Machine.GetPoint( id3 );
    curveLine->start_p = * p1;
    curveLine->end_p = p2;
    curveLine->center_p = p3;
    curveList.push_back( std::move( curveLine ) );
}

void CurveMachine::AddParabolic( int id1, int id2 )
{
    auto curveLine = std::make_unique< CurveLine >();
    curveLine->lineType = PARABOLIC;
    PointType * p1 = point_Machine.GetPoint( id1 );
    PointType * p2 = point_Machine.GetPoint( id2 );
    curveLine->start_p = * p1;
    curveLine->end_p = * p2;
    curveList.push_back( std::move( curveLine ) );
}

CurveLine * CurveMachine::GetCurve( int curveId )
{
    return this->curveList[ curveId ].get();
}

EndNameSpace
