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
#include "NsLusgs.h"
BeginNameSpace( ONEFLOW )

class UNsLusgs : public NsLusgs
{
public:
    UNsLusgs ();
    ~UNsLusgs();
public:
    void LowerSweep();
    void UpperSweep();
    void SingleSweep();
    void PrepareSweep();
    void Update();

    void SolveLowerCell();
    void SolveUpperCell();

    void SolveLower( int fId );
    void SolveUpper( int fId );

    bool CanNotLowerSolve( int fId );
    bool CanNotUpperSolve( int fId );
    void Solve( int fId, int signValue );
    void SetMeshGeometry();
    void PrepareData();
    void PrepareDataFacePrim();
    void CalcViscousTerm();
    void Init();
    void CalcSpectrum();
};

EndNameSpace