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

#include "TurbUpdate.h"
#include "TurbSaPolicy.h"
#include "TurbSstPolicy.h"
#include <memory>
#include "UTurbUpdate.h"
#include "NsInvFlux.h"
#include "TurbCom.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

TurbUpdate::TurbUpdate()
{
}

TurbUpdate::~TurbUpdate()
{
}

void TurbUpdate::CalcFlowField2Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.q[ iEqu ] = SstDensityScaledUpdate(
            turbcom.q0[ iEqu ], turbcom.dq[ iEqu ], turbcom.rho );
    }

    turbcom.q[ IKE ] = SstFloorKineticEnergy( turbcom.q[ IKE ], turbcom.kelim );
    turbcom.q[ IKW ] = SstFloorOmega( turbcom.q[ IKW ], turbcom.kwlim );

    if ( turbcom.transition_model == ITReGama )
    {
        turbcom.q[ ITGama ] = SstCapIntermittency(
            turbcom.q[ ITGama ], turbcom.kGamaLim );
        turbcom.q[ ITRect ] = SstFloorReTheta(
            turbcom.q[ ITRect ], turbcom.kRectLim );
    }
}

void TurbUpdate::CalcFlowField1Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.q[ iEqu ] = SaUnderRelaxedUpdate(
            turbcom.q0[ iEqu ], turbcom.dq[ iEqu ], turbcom.crelax );
    }
}

void TurbUpdate::ModifyValue1Equ()
{
    if ( turbcom.q[ ISA ] < 0.0 )
    {
        ++ turbcom.nneg;

        this->SmoothTurbulencePoint();

        turbcom.q[ ISA ] = SaNonNegativeWorkingVariable(
            turbcom.q[ ISA ], turbcom.inflow[ ISA ] );
    }

    turbcom.xsi   = SaChi( ABS( turbcom.q[ ISA ] ), turbcom.rho / ( turbcom.visl + SMALL ) );
    turbcom.xsi3  = POWER3( turbcom.xsi );
    turbcom.fv1  = SaFv1( turbcom.xsi, turbcom.cv13 );

    const Real nuMax = SaWorkingVariableUpperBound(
        turbcom.maxvistratio, turbcom.visl, turbcom.rho, turbcom.fv1 );

    if ( turbcom.q[ ISA ] > nuMax )
    {
        this->SmoothTurbulencePoint();

        if ( turbcom.q[ ISA ] > nuMax )
        {
            turbcom.q[ ISA ] = nuMax;
        }

        ++ turbcom.npos;
    }

    if ( turbcom.maxvist < turbcom.q[ ISA ] )
    {
        turbcom.maxvist = turbcom.q[ ISA ];
    }
}

std::unique_ptr<Update> CreateTurbUpdate()
{
    return std::make_unique<UTurbUpdate>();
}

EndNameSpace
