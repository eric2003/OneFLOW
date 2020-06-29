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

#include "Stress.h"
#include "Force.h"
#include "Constant.h"

BeginNameSpace( ONEFLOW )

Stress stress;

Stress::Stress()
{
    ;
}

Stress::~Stress()
{
    ;
}

void Stress::CalcStress()
{
    Real divv2p3 = two3rd * ( dudx + dvdy + dwdz );

    txx = viscosity * ( two * dudx - divv2p3 );
    tyy = viscosity * ( two * dvdy - divv2p3 );
    tzz = viscosity * ( two * dwdz - divv2p3 );
    txy = viscosity * ( dudy + dvdx );
    txz = viscosity * ( dudz + dwdx );
    tyz = viscosity * ( dvdz + dwdy );
}

void Stress::CalcForce( Force * force )
{
    this->CalcStress();

    Real coef = - two * area * orey;

    force->x = coef * ( fnx * txx + fny * txy + fnz * txz );
    force->y = coef * ( fnx * txy + fny * tyy + fnz * tyz );
    force->z = coef * ( fnx * txz + fny * tyz + fnz * tzz );
}

EndNameSpace