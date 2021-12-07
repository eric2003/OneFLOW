/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

class ScalarZoneTmp
{
public:
    ScalarZoneTmp();
    ~ScalarZoneTmp();
public:
    void Init( int ist, int ied );
    void InitField( std::vector< double > & uGlobal );
    void Solve( double coef );
    void UpdateUN();
    void GatherField( std::vector< double > & ugfield );
    void CompareField( std::vector< double > & uGlobal );
    void SetBc( int bcL, int bcR );
    double GetRightBcValue();
    void SetLeftBcValue( double lv );
public:
    int zoneid;
    int ist, ied;
    int nNodes;
    std::vector< double > u;
    std::vector< double > un;
    std::vector< double > x;

    std::vector< int > bcPointList;
    std::vector< int > bcTypeList;
};

class ScalarPara;
class ScalarZoneTmp;
class ScalarGrid;

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
    void SolvePart( ScalarZoneTmp * scalarZone );
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

    std::vector< double > u;
    std::vector< double > un;
    std::vector< double > x;
    std::vector< ScalarZoneTmp * > scalarZones;
    ScalarGrid * scalarGrid;
public:
    std::vector< double > du, dua;
    std::vector< double > utheory;
    double l1Norm, l2Norm;
};

EndNameSpace
