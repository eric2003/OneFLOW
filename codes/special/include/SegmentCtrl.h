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
#include "PointMachine.h"

BeginNameSpace( ONEFLOW )

class SegmentCopy;
class FileIO;

class SegmentCtrl
{
public:
    SegmentCtrl();
    ~SegmentCtrl();
public:
    int id;
    int nPoint;
    Real ds1, ds2, lenth;
    int distribution;
    int c1, c2;
    SegmentCopy * segmentCopy;
public:
    Real cA1, cA2, cA3, cA4;
    Real expCoeff, ratio1, ratio2, ratio;
public:
    RealField factorList;
    IntField pidxList;
public:
    void Read( FileIO * ioFile );
    void CalcFactor();
    Real CalcFactor( Real compCoor );
    Real CalcDFactor( Real compCoor );
    Real CalcSFactor( Real compCoor );
    void InitDExp();
    void InitCoef();
    void Init();
    Real CalDExp();
    void SetPara( Real diff, Real & v_min, Real & v_max );
    void CalcEffectiveRatioTest();
    void CalcEffectiveRatio();
    void CalcEffectiveFactor();
    void CopyFactor();
};

class SegmentCopy
{
public:
    SegmentCopy();
    ~SegmentCopy();
public:
    IntField lineList;
    IntField lowerIdList, upperIdList;
public:
    void Read( FileIO * ioFile );
    int GetNSegment() { return lineList.size(); }
    int GetNPoint();
    PointType & GetPoint( int id, int signFlag = 1 );
    void FindSegmentAndLocalId( int id, int signFlag, int & segmentId, int & localId );
};

void GetExponentialCoeff( Real targetRatio, Real cc, Real & coef );

EndNameSpace