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

class TurbLusgsData
{
public:
    TurbLusgsData();
    ~TurbLusgsData();
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

    RealField matrix;
    RealField radius;

    RealField dqj;    //dq of neighbor cell
    RealField dqi;    //dq of this unit
    RealField dqi0;   //The previous DQ of this unit (used for normal)
    RealField primj;   //�ھӵ�Ԫ��ԭʼ����
    RealField primF;   //��Ԫ�ϵ�ԭʼ����
    RealField rhs0;   //��ʽ�в�����
    RealField dfj;   //��ʱ���飬��������ʽ�в�����
    RealField drhs;   //����nsweep>1�����
    RealField rhs ;   //����nʱ�̵��Ҷ���
    RealField tmp; //��ʱ����
};

extern TurbLusgsData turblu;

class TurbLusgs : public LusgsSolver
{
public:
    TurbLusgs ();
    ~TurbLusgs();
public:
    void InitializeSub();
public:
    void DumpSweepInformation();
    void ZeroFluxIncrement   ();
    void AddFluxIncrement    ();
    void AddFluxIncrement( const Real & coef );
    void GetFluxIncrement( int signOfMatrix );
    void GetStandardFluxIncrement( int signOfMatrix );
    void InitializeSweep( int iSweep );
    bool UpdateSweep    ( int iSweep );
public:
    void CmpLowerChange();
    void CmpUpperChange();
    bool IsOversetCell  ();
    void ZeroOversetCell();
};


EndNameSpace