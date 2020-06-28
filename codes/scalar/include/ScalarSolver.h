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
#include "Configure.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

const int SCALAR_INTERFACE   = -1;
const int SCALAR_COMPUTE     = 1;
const int SCALAR_EXTRAPOLATE = 2;

class ScalarZone
{
public:
    ScalarZone();
    ~ScalarZone();
public:
    void Init( int ist, int ied );
    void InitField( vector< double > & uGlobal );
    void Solve( double coef );
    void UpdateUN();
    void GatherField( vector< double > & ugfield );
    void CompareField( vector< double > & uGlobal );
    void SetBc( int bcL, int bcR );
    double GetRightBcValue();
    void SetLeftBcValue( double lv );
public:
    int zoneid;
    int ist, ied;
    int nNode;
    vector< double > u;
    vector< double > un;
    vector< double > x;

    vector< int > bcPointList;
    vector< int > bcTypeList;
};

class ScalarPara;
class ScalarZone;
class ScalarSolver
{
public:
    ScalarSolver();
    ~ScalarSolver();
public:
    void Run();
    void RunTest( ScalarPara * para );
    void Init();
    void Solve();
    void SolvePart( ScalarZone * scalarZone );
    void SolveFlowField();
    void Boundary();
    void InitCtrlParameter();
    void InitCtrlParameterTest( ScalarPara * para );
    void InitGrid();
    void InitFlowField();
    void InitZoneFlowField();
    void UpdateUN();
    void CompareField();
public:
    void SetScalarZone();
    void SolveOneStep();
    void Visual();
    void Visual( ScalarPara * para );
    void Theory( double t );
    double ScalarFun( double xm );
    double SquareFun( double xm );
    double SinFun( double xm );
    double MyCurve( double xm );
    void FreeScalarZones();
public:
    void CalcL1Norm();
    void CalcL2Norm();
public:
    int nx;

    double len;
    double dx;
    int    nt;
    double dt;
    double c;

    vector< double > u;
    vector< double > un;
    vector< double > x;
    vector< ScalarZone * > scalarZones;
public:
    vector< double > du, dua;
    vector< double > utheory;
    double l1Norm, l2Norm;
};

EndNameSpace