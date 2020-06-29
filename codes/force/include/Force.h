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

BeginNameSpace( ONEFLOW )

class Force
{
public:
    Force();
    ~Force();
    Force( const Force & rhs );
    Force( const Real & value );
public:
    Real x, y, z;
public:
    Force & operator += ( const Force & rhs );
    Force & operator = ( const Real & value );
    Force & operator = ( const Force & rhs );
    Force & operator /= ( const Real & value );
};

Force operator + ( const Force & f1, const Force & f2 );
Force operator / ( const Force & f, const Real & value );

class AeroCom
{
public:
    AeroCom();
    ~AeroCom();
public:
    void Init();
    Real CalcCL( Force * f );
    Real CalcCD( Force * f );
    Real CalcCF( Force * f, Real area );
public:
    Real aref, lref;
    Real xref, yref, zref;
    Real aoa, aos;
    Real sina, cosa;
    Real sinb, cosb;
    Real cForce;
    Real cMoment;
};


class AeroForce
{
public:
    AeroForce ();
    ~AeroForce();
public:
    Force aero;
    Force pres, vis;
    Force total;
    Force aeromom, mom;

    Real cl, cd;
    Real cp, cf, cdCl2Pa;
    Real vfx, vfy, vfz;
    Real power;
public:
    void SumForce();
    void CalcPower();
    void CalcMoment( Real xc, Real yc, Real zc );
    void AddForce( AeroForce * rhs );
    void Init();
};

class AeroForceInfo
{
public:
    AeroForceInfo();
    ~AeroForceInfo();
public:
    void Init();
    void CollectForce();
    void Sum( Force * force );
    void Sum( Real * var );
    void CalcCoef();
public:
    AeroForce totalForce;
    Force cf, cpres, cmom;
    Real cpower;
    Real cd, cl, cd_pres, cd_vis;
    Real cdl;
    Real pres_center;
};

extern AeroForceInfo aeroForceInfo;
extern AeroCom aeroCom;

EndNameSpace