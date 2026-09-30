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
#include "HXLookup.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class SegmentCtrl;
class CurveInfo;
class CurveMesh;
class TextFileParser;

class LineMachine
{
public:
    LineMachine();
    ~LineMachine();
public:
    HXVector< std::unique_ptr< SegmentCtrl > > segmentCtrlList;
    HXVector< std::unique_ptr< CurveInfo > > curveInfoList;
    HXVector< CurveMesh * > curveMeshList;
    IntField dimList;
    RealField ds1List, ds2List;
public:
    HXLookup<int> lineLookup;
    LinkField lineList; 
public:
    int AddLine( int p1, int p2 );
    void AddLine( int p1, int p2, int id );
    void AddCircle( int p1, int pc, int p2, int id );
    void AddDimension( TextFileParser & textFileParser );
    void AddDs( TextFileParser & textFileParser );
    void GenerateAllLineMesh();
    void CreateAllLineMesh();
public:
    SegmentCtrl * GetSegmentCtrl( int id ) const;
    CurveMesh * GetCurveMesh( int id ) const;
    CurveInfo * GetCurveInfo( int id ) const;
public:
    CurveMesh * GetLineMeshByTwoPoint( const int & p1, const int & p2, int & direction ) const;
    int GetLineIdByTwoPoint( const int & p1, const int & p2 ) const;
};

extern LineMachine line_Machine;

EndNameSpace
