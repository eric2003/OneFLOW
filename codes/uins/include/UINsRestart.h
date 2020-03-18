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
#include "Restart.h"
BeginNameSpace( ONEFLOW )

class UINsRestart : public Restart
{
public:
    UINsRestart();
    ~UINsRestart();
public:
    void InitinsRestart( int sTid );
};

class ShockVertex
{
public:
    ShockVertex();
    ~ShockVertex();
public:
    Real xc, yc;
    Real Ms, gama;
    Real ru, uu, vu, pu;
    Real rd, ud, vd, pd;
public:
    void Init();
    Real vortexfun( Real x, Real y, Real gama );
    void Cal( Real x, Real y, Real gama, Real & vx, Real & vy );
};

EndNameSpace