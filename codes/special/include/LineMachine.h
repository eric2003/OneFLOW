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
#include "Mid.h"

BeginNameSpace( ONEFLOW )

class SegmentCtrl;
class CurveInfo;
class CurveMesh;
class FileIO;

class LineMachine
{
public:
    LineMachine();
    ~LineMachine();
public:
    HXVector< SegmentCtrl * > segmentCtrlList;
    HXVector< CurveInfo * > curveInfoList;
    HXVector< CurveMesh * > curveMeshList;
    IntField dimList;
    RealField ds1List, ds2List;
public:
    set< Mid<int> > refLines;
    LinkField lineList; 
public:
    int AddLine( int p1, int p2 );
    void AddLine( int p1, int p2, int id );
    void AddDimension( FileIO * ioFile );
    void AddDs( FileIO * ioFile );
    void GenerateAllLineMesh();
    void CreateAllLineMesh();
public:
    SegmentCtrl * GetSegmentCtrl( int id );
    CurveMesh * GetCurveMesh( int id );
    CurveInfo * GetCurveInfo( int id );
public:
    CurveMesh * GetLineMeshByTwoPoint( const int & p1, const int & p2, int & direction );
    int GetLineIdByTwoPoint( const int & p1, const int & p2 );
};

extern LineMachine line_Machine;

EndNameSpace