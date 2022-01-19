/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "HXType.h"
#include "ScalarGrid.h"
#include <vector>


BeginNameSpace( ONEFLOW )

class FieldPara;

class ScalarField
{
public:
    ScalarField();
    ~ScalarField();
public:
    ScalarGrid * grid;
    RealList q;
    RealList res;
    RealList qL, qR;
    RealList invflux;
    FieldPara * para;
    Real qInf;
public:
    void InitFlowField( ScalarGrid * grid );
    Real ScalarFun( Real xm );
    Real SquareFun( Real xm );
public:
    void ToTecplot( RealList & varlist, std::string const & fileName );
    void SolveFlowField( FieldPara * para );
    void UpdateResidual();
    void TimeIntergral();
    void SolveOneStep();
    void Update();
    void Boundary();
    void GetQLQR();
    void CalcInvFlux();
    void AddF2CField( RealList & cField, RealList & fField );
    void SetParaPointer( FieldPara * para );
public:
    void Visual();
    void AddVisualData( RealList & qList, RealList & theoryList, RealList & xcoorList );
    void Theory( Real time, RealList & theory );
};


EndNameSpace
