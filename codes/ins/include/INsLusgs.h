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
#include "Lusgs.h"
BeginNameSpace( ONEFLOW )

class ILusgsData
{
public:
    ILusgsData();
    ~ILusgsData();
public:
    void Init();
public:
    int  nEqu;
    int  nBEqu;
    int  numberOfSweeps;
    int  numberOfRealSweeps;
    bool keyPrim;
    Real tol;
    Real norm0, dqSweep, dmax;
    Real norm;
    Real gama, rad;

    Real lmdOnFace1;
    Real lmdOnFace2;
    Real lmdOnFace3;

    Real visrad;

    RealField radius;
    RealField dqj;    //dq of neighbor cell
    RealField dqi;    //dq of this unit
    RealField dqi0;   //The previous DQ of this unit (used for normal)
    RealField primj;   //邻居单元的原始变量
    RealField primF;   //面元上的原始变量
    RealField rhs0;   //隐式残差增量
    RealField dfj;   //临时数组，用于求隐式残差增量
    RealField drhs;   //用于nsweep>1的情况
    RealField rhs ;   //方程n时刻的右端项
    RealField tmp; //临时数组
};

extern ILusgsData nslu;

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
    void CmpFaceEigenValue( RealField & prim );
    void GetStandardFluxIncrement( int signOfMatrix );
    void InitializeSweep( int iSweep );
    bool UpdateSweep    ( int iSweep );
public:
    void CmpLowerChange();
    void CmpUpperChange();
    bool IsOversetCell  ();
    void ZeroOversetCell();
};

void CmpIDH( RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy );


EndNameSpace