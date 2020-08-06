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

#include "Sutherland.h"
#include <cmath>
using namespace std;

BeginNameSpace( ONEFLOW )

Sutherland::Sutherland()
{
    this->SetInvalidValue();
}

Sutherland::~Sutherland()
{
}

void Sutherland::SetInvalidValue()
{
    this->t0_dim = -1.0;
    this->ts_dim = -1.0;
    this->mu0_dim = -1.0;

    this->ts = -1.0;

}

void Sutherland::Init( Real tref_dim )
{
    this->t0_dim = 273.16; // zero degree celsius temperature(K)
    this->ts_dim = 110.4; //dimensional temperature constant in sutherland formula
    this->mu0_dim = 1.715e-5; //dimensional air viscosity at zero celsius degree

    this->ts = this->ts_dim / tref_dim;
}

Real Sutherland::CalcViscosityDim( Real t_dim )
{
    Real coef = ( t0_dim + ts_dim ) / ( t_dim + ts_dim );
    Real vis_dim = coef * pow( t_dim / t0_dim, 1.5 ) * this->mu0_dim;
    return vis_dim;
}

Real Sutherland::CalcViscosity( Real t )
{
    Real t3 = t * t * t;
    return sqrt( t3 ) * ( 1 + this->ts ) / ( t + this->ts );
}

EndNameSpace