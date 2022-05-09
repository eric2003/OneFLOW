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
#include "Lusgs.h"
BeginNameSpace( ONEFLOW )

class LusgsData
{
public:
    LusgsData();
    ~LusgsData();
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
    RealField primj;   //The original variable of neighbor unit
    RealField primF;   //The original variable on the surface
    RealField rhs0;   //Implicit residual increment
    RealField dfj;   //Temporary array,Used to findImplicit residual increment
    RealField drhs;   //For nsweep > 1
    RealField rhs ;   //Of the equation at time nRight end item
    RealField tmp; //Temporary array
};

extern LusgsData nslu;

class NsLusgs : public LusgsSolver
{
public:
    NsLusgs ();
    ~NsLusgs();
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

void CalcDH( RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy );


EndNameSpace
