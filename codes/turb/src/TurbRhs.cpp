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

#include "TurbRhs.h"
#include "Ctrl.h"
#include "UTurbInvFlux.h"
#include "UTurbVisFlux.h"
#include "UTurbSrcFlux.h"
#include "UTurbSpectrum.h"
#include "UTurbUnsteady.h"
#include "UTurbBcSolver.h"

BeginNameSpace( ONEFLOW )

TurbRhs::TurbRhs()
{
    ;
}

TurbRhs::~TurbRhs()
{
    ;
}

void TurbRhs::CmpRHS()
{
    TurbCmpRHS();
}

void TurbCmpBc()
{
    UTurbBcSolver * uTurbBcSolver = new UTurbBcSolver();
    uTurbBcSolver->Init();
    uTurbBcSolver->CmpBc();
    delete uTurbBcSolver;
}

void TurbCmpRHS()
{
    TurbCmpBc();

    TurbCmpSpectrum();

    TurbCmpSrcFlux();

    TurbCmpInvFlux();

    TurbCmpVisFlux();

    TurbCmpDualTimeStepSrc();
}


void TurbCmpInvFlux()
{
    UTurbInvFlux * uTurbInvFlux = new UTurbInvFlux();
    uTurbInvFlux->CmpFlux();
    delete uTurbInvFlux;
}

void TurbCmpVisFlux()
{
    UTurbVisFlux * uTurbVisFlux = new UTurbVisFlux();
    uTurbVisFlux->CmpVisFlux();
    delete uTurbVisFlux;
}

void TurbCmpSrcFlux()
{
    UTurbSrcFlux * uTurbSrcFlux = new UTurbSrcFlux();
    uTurbSrcFlux->CmpSrcFlux();
    delete uTurbSrcFlux;
}

void TurbCmpSpectrum()
{
    UTurbSpectrum * uTurbSpectrum = new UTurbSpectrum();
    uTurbSpectrum->CmpSpectrum();
    delete uTurbSpectrum;
}

void TurbCmpDualTimeStepSrc()
{
    //dual time step source
    if ( ctrl.idualtime == 1 )
    {
        UTurbUnsteady * uTurbUnsteady = new UTurbUnsteady();
        uTurbUnsteady->CmpDualTimeSrc();
        delete uTurbUnsteady;
    }
}

void CmpTurbulentViscosity()
{
    UTurbSrcFlux * uTurbSrcFlux = new UTurbSrcFlux();
    uTurbSrcFlux->CmpVist();
    delete uTurbSrcFlux;

}

EndNameSpace