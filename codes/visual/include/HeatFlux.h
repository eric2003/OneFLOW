/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include <memory> // Added for std::unique_ptr

BeginNameSpace( ONEFLOW )

class SurfaceValue
{
public:
    SurfaceValue();
    ~SurfaceValue();

public:
    // FIX: Replaced raw pointer with std::unique_ptr for automatic memory management.
    std::unique_ptr<RealField> var;
};

class HeatFlux
{
public:
    HeatFlux ();
    ~HeatFlux();

public:
    // FIX: Replaced raw pointer arrays with std::unique_ptr arrays.
    HXVector< std::unique_ptr<SurfaceValue> > heatflux;
    HXVector< std::unique_ptr<SurfaceValue> > fricflux;
    IntField flag;
    bool init_flag;

public:
    void Init();
    void InitGlobal();
    void Allocate();
    void DeAllocate();
};

extern HeatFlux heat_flux;

EndNameSpace
