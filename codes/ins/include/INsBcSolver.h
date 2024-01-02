/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

class INsBcSolver
{
public:
    INsBcSolver();
    ~INsBcSolver();
public:
    typedef void ( INsBcSolver:: * BcPointer )();
public:
    void SetBc();
    void SetSolidSurfaceBc();
    virtual void Solve(){};
    BcPointer bcPointer;
    bool updateFlag;
public:
    void InFlowBc           ();
    void OutFlowBc          ();
    void FarFieldBc         ();
    void IsothermalVisWallBc();
    void AdiabaticVisWallBc ();
    void SymmetryBc         ();
    void PoleBc             ();
    void OversetBc          ();
    void InterfaceBc        ();
    void NoBc               ();
    void UserDefinedBc      ();
    void PeriodicBc         ();
    void VelocityBc();
public:
    void CalcFaceBc();
    Real CalcReciMolecularWeight( RealField & prim );
    Real CalcDensity( RealField & prim, Real pres, Real temperature );
};

class BcData;
extern BcData ins_bc_data;

EndNameSpace
