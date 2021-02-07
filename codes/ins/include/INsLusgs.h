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
#include "Lusgs.h"
BeginNameSpace( ONEFLOW )

class INsLusgs : public LusgsSolver
{
public:
    INsLusgs ();
    ~INsLusgs();
public:
    void InitializeSub();
public:
    void DumpSweepInformation();
    void ZeroFluxIncrement   ();
    void AddViscousTerm      ();
    void AddFluxIncrement    ();
    void AddFluxIncrement( const Real & coef );
    void GetFluxIncrement( int signOfMatrix );
    void CalcFaceEigenValue( RealField & prim );
    void GetStandardFluxIncrement( int signOfMatrix );
    void InitializeSweep( int iSweep );
    bool UpdateSweep    ( int iSweep );
public:
    void CalcLowerChange();
    void CalcUpperChange();
    bool IsOversetCell  ();
    void ZeroOversetCell();
};

void CalcIDH( RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy );


EndNameSpace