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

#include "UsdBasic.h"
#include "Ctrl.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UsdBasic::UsdBasic()
{
    ;
}

UsdBasic::~UsdBasic()
{
    ;
}

void UsdBasic::InitCoef()
{
    // linearTwoStepMethods
    // 0  - Euler explicit                          : Order 1
    // 1  - Bakward Euler                           : Order 1   A-stable
    // 2  - One-step trapezoidal( Crank-Nicolson )  : Order 2   A-stable
    // 3  - Bakward differentiation                 : Order 2   A-stable
    // 4  - Adams type                              : Order 2   A-stable
    // 5  - Lees type                               : Order 2   A-stable
    // 6  - Two-step trapezoidal                    : Order 2   A-stable
    // 7  - Leapfrog                                : Order 2
    // 8  - Adams-Bashforth                         : Order 2
    // 9  - Third-order implicit                    : Order 3
    // 10 - Adams-Moulton                           : Order 3
    // 11 - Milne                                   : Order 4

    //     theat   ksai    phi   name
    // 0   0       0       0
    // 1   1       0       0
    // 2   1/2     0       0
    // 3   1       1/2     0
    // 4   3/4     0      -1/4
    // 5   1/3    -1/2    -1/3
    // 6   1/2    -1/2    -1/2
    // 7   0      -1/2     0
    // 8   0       0      1/2
    // 9   1/3    -1/6     0
    // 10  5/12    0      1/12
    // 11  1/ 6   -1/2    -1/6

    this->coeff.resize( 3 );

    Real thet = zero;
    Real xi   = zero;
    Real phi  = zero;

    int linearTwoStepMethods = ctrl.linearTwoStepMethods;

    if ( linearTwoStepMethods == 0 )
    {
        // 0  - Euler explicit                          : Order 1
        // 0   0       0       0
        thet = zero;
        xi   = zero;
        phi  = zero;
    }
    else if ( linearTwoStepMethods == 1 )
    {
        // 1  - Bakward Euler                          : Order 1   A-stable
        // 1   1       0       0
        thet = 1.0;
        xi   = zero;
        phi  = zero;
    }
    else if ( linearTwoStepMethods == 2 )
    {
        // 2  - One-step trapezoidal( Crank-Nicolson )  : Order 2   A-stable
        // 2   1/2     0       0
        thet = 0.5;
        xi   = zero;
        phi  = zero;
    }
    else if ( linearTwoStepMethods == 3 )
    {
        // 3  - Bakward differentiation                : Order 2   A-stable
        // 3   1       1/2     0
        thet = 1.0;
        xi   = 0.5;
        phi  = zero;
    }
    else if ( linearTwoStepMethods == 4 )
    {
        // 4  - Adams type                              : Order 2   A-stable
        // 4   3/4     0      -1/4
        thet = 3.0 / 4.0;
        xi   = 0.0;
        phi  = - 1.0 / 4.0;
    }
    else if ( linearTwoStepMethods == 5 )
    {
        // 5  - Lees type                               : Order 2   A-stable
        // 5   1/3    -1/2    -1/3
        thet =   1.0 / 3.0;
        xi   = - 1.0 / 2.0;
        phi  = - 1.0 / 3.0;
    }
    else if ( linearTwoStepMethods == 6 )
    {
        // 6  - Two-step trapezoidal                    : Order 2   A-stable
        // 6   1/2    -1/2    -1/2
        thet =   1.0 / 2.0;
        xi   = - 1.0 / 2.0;
        phi  = - 1.0 / 2.0;
    }
    else if ( linearTwoStepMethods == 7 )
    {
        // 7  - Leapfrog                                : Order 2
        // 7   0      -1/2     0
        thet =   0.0;
        xi   = - 1.0 / 2.0;
        phi  =   0.0;
    }
    else if ( linearTwoStepMethods == 8 )
    {
        // 8  - Adams-Bashforth                         : Order 2
        // 8   0       0      1/2
        thet =   0.0;
        xi   =   0.0;
        phi  =   1.0 / 2.0;
    }
    else if ( linearTwoStepMethods == 9 )
    {
        // 9  - Third-order implicit                    : Order 3
        // 9   1/3    -1/6     0
        thet =   1.0 / 3.0;
        xi   = - 1.0 / 6.0;
        phi  =   0.0;
    }
    else if ( linearTwoStepMethods == 10 )
    {
        // 10 - Adams-Moulton                           : Order 3
        // 10  5/12    0      1/12
        thet =   5.0 / 12.0;
        xi   =   0.0;
        phi  =   1.0 / 12.0;
    }
    else if ( linearTwoStepMethods == 11 )
    {
        // 11 - Milne                                   : Order 4
        // 11  1/ 6   -1/2    -1/6
        thet =   1.0 / 6.0;
        xi   = - 1.0 / 2.0;
        phi  = - 1.0 / 6.0;
    }
    else
    {
        cout << " Error !!!!! linearTwoStepMethods = " << linearTwoStepMethods << endl;
    }

    coeff[ 0 ] = thet;
    coeff[ 1 ] = xi  ;
    coeff[ 2 ] = phi ;
}

void UsdBasic::CalcResCoef()
{
    Real thet = coeff[ 0 ];
    Real xi   = coeff[ 1 ];
    Real phi  = coeff[ 2 ];

    if ( ABS( thet ) < 1.0e-3 )
    {
        resc1 = thet;
        resc2 = 1.0 - thet + phi;
        resc3 = - phi;
    }
    else
    {
        resc1 = 1.0;
        resc2 = ( 1.0 - thet + phi ) / thet;
        resc3 = - phi / thet;
    }
}

void UsdBasic::CalcSpectrumCoeff()
{
    Real thet = coeff[ 0 ];
    Real xi   = coeff[ 1 ];
    Real phi  = coeff[ 2 ];

    Real coef;

    if ( ABS( thet ) < 1.0e-3 )
    {
        coef = 1.0;
    }
    else
    {
        coef = 1.0 / ( thet );
    }

    sp1 = coef * 1.0;
    sp2 = coef * ( 1.0 + xi );
}

void UsdBasic::CalcSrcCoeffBasic()
{
    Real thet = coeff[ 0 ];
    Real xi   = coeff[ 1 ];
    Real phi  = coeff[ 2 ];
    
    if ( ABS( thet ) < 1.0e-3 )
    {
        bsc1 = ( 1.0 + xi );
        bsc3 = ( xi       );
    }
    else
    {
        bsc1 = ( 1.0 + xi ) / thet;
        bsc3 = ( xi       ) / thet;
    }
    bsc2 = - bsc1 - bsc3;
}

void UsdBasic::CalcSrcCoeff()
{
    this->CalcSrcCoeffBasic();

    sc1 = bsc1 / ctrl.pdt;
    sc3 = bsc3 / ctrl.pdt1;
    sc2 = - sc1 - sc3;
}

void UsdBasic::InitBasic()
{
    this->InitCoef();
    this->CalcResCoef();
    this->CalcSpectrumCoeff();
}

EndNameSpace