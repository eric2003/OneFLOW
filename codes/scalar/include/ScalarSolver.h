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

class ScalarZone
{
public:
    ScalarZone();
    ~ScalarZone();
public:
    int ist, ied;
};

class ScalarPara;

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
    void SolvePart();
    void SolvePart( int ist, int ied );
    void SolveFlowField();
    void InitCtrlParameter();
    void InitCtrlParameterTest( ScalarPara * para );
    void InitGrid();
    void InitFlowField();
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
    void ComputeL1Norm();
    void ComputeL2Norm();
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
    vector< double > u0;
    vector< double > du, dua;
    vector< double > utheory;
    double l1Norm, l2Norm;
};

EndNameSpace