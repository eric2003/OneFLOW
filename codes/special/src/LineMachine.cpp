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

#include "LineMachine.h"
#include "SegmentCtrl.h"
#include "CurveInfo.h"
#include "LineInfo.h"
#include "LineMesh.h"
#include "LineMeshImp.h"
#include "FileIO.h"
#include "HXMath.h"
#include <iostream>
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

LineMachine line_Machine;

LineMachine::LineMachine()
{
}

LineMachine::~LineMachine()
{
    for ( int i = 0; i < curveInfoList.size(); ++ i )
    {
        delete curveInfoList[ i ];
    }

    for ( int i = 0; i < segmentCtrlList.size(); ++ i )
    {
        delete segmentCtrlList[ i ];
    }
}

SegmentCtrl * LineMachine::GetSegmentCtrl( int id )
{
    int idx = ABS( id ) - 1;
    return this->segmentCtrlList[ idx ];
}

CurveMesh * LineMachine::GetCurveMesh( int id )
{
    int idx = ABS( id ) - 1;
    return this->curveMeshList[ idx ];
}

CurveInfo * LineMachine::GetCurveInfo( int id )
{
    int idx = ABS( id ) - 1;
    return this->curveInfoList[ idx ];
}

int LineMachine::AddLine( int p1, int p2 )
{
    IntField line;
    line.push_back( p1 );
    line.push_back( p2 );

    int n = this->refLines.size();

    Mid<int> fMid( 2, n + 1 );
    fMid.data = line;
    std::sort( fMid.data.begin(), fMid.data.end() );

    set< Mid<int> >::iterator iter = this->refLines.find( fMid );
    if ( iter == this->refLines.end() )
    {
        this->refLines.insert( fMid );
        this->lineList.push_back( line );
        return n + 1;
    }
    else
    {
        return iter->id;
    }

}

void LineMachine::AddLine( int p1, int p2, int id )
{
    int idd = this->AddLine( p1, p2 );
    CurveInfo * line = new LineInfo( p1, p2, id );
    this->curveInfoList.push_back( line );

    SegmentCtrl * segmentCtrl = new SegmentCtrl();
    segmentCtrl->id = id;
    this->segmentCtrlList.push_back( segmentCtrl );
}

void LineMachine::AddDimension( FileIO * ioFile )
{
    int id = ioFile->ReadNextDigit< int >();
    int dim = ioFile->ReadNextDigit< int >();
    this->dimList.push_back( dim );
    SegmentCtrl * segmentCtrl = this->GetSegmentCtrl( id );
    segmentCtrl->nPoint = dim;
}

void LineMachine::AddDs( FileIO * ioFile )
{
    int id = ioFile->ReadNextDigit< int >();
    SegmentCtrl * segmentCtrl = this->GetSegmentCtrl( id );
    segmentCtrl->Read( ioFile );
}

void LineMachine::CreateAllLineMesh()
{
    int nLine = curveInfoList.size();
    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ];
        int lineType = curveInfo->type;
        CurveMesh * curveMesh = CreateLineMesh( curveInfo );
        curveMesh->segmentCtrl = this->GetSegmentCtrl( curveInfo->id );
        this->curveMeshList.push_back( curveMesh );
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
            CurveMesh * curveMesh = this->curveMeshList[ iLine ];
            curveMesh->GenerateLineMesh();

            if ( curveMesh->state == 1 ) nCount ++;
        }
        if ( nCount == nLine ) break;
    }
}

CurveMesh * LineMachine::GetLineMeshByTwoPoint( const int & p1, const int & p2, int & direction )
{
    direction = 1;
    int nLine = curveInfoList.size();

    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ];
        if ( curveInfo->p1 == p1 &&
             curveInfo->p2 == p2 )
        {
            direction = 1;
            int & lineId = curveInfo->id;
            return this->GetCurveMesh( lineId );
        }
        else if ( curveInfo->p2 == p1 &&
                  curveInfo->p1 == p2 )
        {
            direction = - 1;
            int & lineId = curveInfo->id;
            return this->GetCurveMesh( lineId );
        }
    }
    return 0;
}

int LineMachine::GetLineIdByTwoPoint( const int & p1, const int & p2 )
{
    int nLine = curveInfoList.size();

    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        CurveInfo * curveInfo = curveInfoList[ iLine ];
        if ( curveInfo->p1 == p1 &&
            curveInfo->p2 == p2 )
        {
            int & lineId = curveInfo->id;
            return lineId;
        }
        else if ( curveInfo->p2 == p1 &&
            curveInfo->p1 == p2 )
        {
            int & lineId = curveInfo->id;
            return lineId;
        }
    }
    return 0;
}


EndNameSpace