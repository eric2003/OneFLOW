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

#include "LineMachine.h"
#include "SegmentCtrl.h"
#include "CurveInfo.h"
#include "LineInfo.h"
#include "CircleInfo.h"
#include "LineMesh.h"
#include "LineMeshImp.h"
#include "TextFileParser.h"
#include "HXMath.h"
#include <iostream>
#include <algorithm>


BeginNameSpace( ONEFLOW )

LineMachine line_Machine;

LineMachine::LineMachine()
{
}

LineMachine::~LineMachine() = default;

SegmentCtrl * LineMachine::GetSegmentCtrl( int id ) const
{
    int idx = ABS( id ) - 1;
    return this->segmentCtrlList[ idx ].get();
}

CurveMesh * LineMachine::GetCurveMesh( int id ) const
{
    int idx = ABS( id ) - 1;
    return this->curveMeshList[ idx ];
}

CurveInfo * LineMachine::GetCurveInfo( int id ) const
{
    int idx = ABS( id ) - 1;
    return this->curveInfoList[ idx ].get();
}

int LineMachine::AddLine(int p1, int p2)
{
    IntField line;
    line.push_back(p1);
    line.push_back(p2);

    //int lineIndex = this->lineLookup.FindOrAdd( line ); // Ensure the line is registered in the lookup
    auto [lineIndex, isNew] = this->lineLookup.FindOrAdd( line ); // Ensure the line is registered in the lookup
    if ( isNew )
    {
        this->lineList.push_back(line);     // Preserve original order
    }
    return lineIndex;
}

void LineMachine::AddLine( int p1, int p2, int id )
{
    this->AddLine( p1, p2 );
    auto line = std::make_unique< LineInfo >( p1, p2, id );
    this->curveInfoList.push_back( std::move( line ) );

    auto segmentCtrl = std::make_unique< SegmentCtrl >();
    segmentCtrl->id = id;
    this->segmentCtrlList.push_back( std::move( segmentCtrl ) );
}

void LineMachine::AddCircle( int p1, int pc, int p2, int id )
{
    this->AddLine( p1, p2 );
    auto circle = std::make_unique< CircleInfo >( p1, pc, p2, id );
    this->curveInfoList.push_back( std::move( circle ) );

    auto segmentCtrl = std::make_unique< SegmentCtrl >();
    segmentCtrl->id = id;
    this->segmentCtrlList.push_back( std::move( segmentCtrl ) );
}

void LineMachine::AddDimension( TextFileParser & textFileParser )
{
    int id = textFileParser.ReadNextDigit< int >();
    int dim = textFileParser.ReadNextDigit< int >();
    this->dimList.push_back( dim );
    SegmentCtrl * segmentCtrl = this->GetSegmentCtrl( id );
    segmentCtrl->nPoint = dim;
}

void LineMachine::AddDs( TextFileParser & textFileParser )
{
    int id = textFileParser.ReadNextDigit< int >();
    SegmentCtrl * segmentCtrl = this->GetSegmentCtrl( id );
    segmentCtrl->Read( & textFileParser );
}

void LineMachine::CreateAllLineMesh()
{
    int nLine = curveInfoList.size();
    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ].get();
        auto curveMesh = CreateLineMesh( curveInfo );
        curveMesh->segmentCtrl = this->GetSegmentCtrl( curveInfo->id );
        this->curveMeshList.push_back( std::move( curveMesh ) );
    }
}

void LineMachine::GenerateAllLineMesh()
{
    CreateAllLineMesh();

    while ( true )
    {
        int nCount = 0;
        int nLine = curveInfoList.size();
        for ( int iLine = 0; iLine < nLine; ++ iLine )
        {
            CurveMesh * curveMesh = this->curveMeshList[ iLine ].get();
            curveMesh->GenerateLineMesh();

            if ( curveMesh->state == 1 ) nCount ++;
        }
        if ( nCount == nLine ) break;
    }
}

CurveMesh * LineMachine::GetLineMeshByTwoPoint( const int & p1, const int & p2, int & direction ) const
{
    direction = 1;
    int nLine = curveInfoList.size();

    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ].get();
        if ( curveInfo->p1 == p1 &&
             curveInfo->p2 == p2 )
        {
            direction = 1;
            const int lineId = curveInfo->id;
            return this->GetCurveMesh( lineId );
        }
        else if ( curveInfo->p2 == p1 &&
                  curveInfo->p1 == p2 )
        {
            direction = - 1;
            const int lineId = curveInfo->id;
            return this->GetCurveMesh( lineId );
        }
    }
    return 0;
}

int LineMachine::GetLineIdByTwoPoint( const int & p1, const int & p2 ) const
{
    int nLine = curveInfoList.size();

    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ].get();
        if ( curveInfo->p1 == p1 &&
            curveInfo->p2 == p2 )
        {
            const int lineId = curveInfo->id;
            return lineId;
        }
        else if ( curveInfo->p2 == p1 &&
            curveInfo->p1 == p2 )
        {
            const int lineId = curveInfo->id;
            return lineId;
        }
    }
    return 0;
}


EndNameSpace
