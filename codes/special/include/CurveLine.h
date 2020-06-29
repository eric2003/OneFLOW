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
#include "Point.h"

BeginNameSpace( ONEFLOW )

typedef Point< Real > PointType;

class CurveLine;

class StrCurveLoop
{
public:
    StrCurveLoop();
    ~StrCurveLoop();
public:
    int ni, nj;
    IntField curveIdList;
public:
    void AddCurve( int curveId );
    CurveLine * GetCurve( int curveId );
    void SetDimension();
};

class CurveLine
{
public:
    CurveLine();
    ~CurveLine();
public:
    RealField x;
    RealField y;
    RealField z;
    int nNode;
public:
    PointType start_p, end_p;
    PointType center_p;
    Real totalCurveLength;
    Real ds1, ds2;
    int lineType;
public:
    void Alloc( int nNode );
    void MakeLine( PointType & p1, PointType & p2 );
    void MakeCircle();
    void MakeCircle( PointType & p1, PointType & p2, PointType & p0 );
    void GenerateCurveLine();
    void GenerateCircleLine();
    void GenerateParabolicLine();
    void CalcTotalLength();
    Real CalcParabolicLength(Real x, Real p);
    void FindXByLength(Real & x, Real p, Real maxX, Real length);
public:
    void CalcLineLength();
    void CalcNormal( RealField & nbx, RealField & nby, RealField & nbz );
};

Real CalcDist( PointType & p1, PointType & p2 );

void CrossProduct( RealField & a, RealField & b, RealField & c );


EndNameSpace