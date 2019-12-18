/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "PointMachine.h"

BeginNameSpace( ONEFLOW )

class CurveInfo;
class SegmentCtrl;

class CurveMesh
{
public:
    CurveMesh();
    virtual ~CurveMesh();
public:
    CurveInfo * curveInfo;
    SegmentCtrl * segmentCtrl;
    HXVector< PointType * > ptList;
    int state;
public:
    bool IsValidState();
    void GenerateCurveMesh();
    PointType & GetPoint( int id, int signFlag );
public:
    int GetDim() { return ptList.size(); };
    virtual void GenerateLineMesh() {};
    virtual void ComputeCurveGeometry() {};
    virtual void ComputeCoor( Real s, Real & xt, Real & yt, Real & zt ) {};
};

void Alloc( RealField2D & field, int ni, int nj );

EndNameSpace