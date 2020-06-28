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

void TurbRhs::CalcRHS()
{
    TurbCalcRHS();
}

void TurbCalcBc()
{
    UTurbBcSolver * uTurbBcSolver = new UTurbBcSolver();
    uTurbBcSolver->Init();
    uTurbBcSolver->CalcBc();
    delete uTurbBcSolver;
}

void TurbCalcRHS()
{
    TurbCalcBc();

    TurbCalcSpectrum();

    TurbCalcSrcFlux();

    TurbCalcInvFlux();

    TurbCalcVisFlux();

    TurbCalcDualTimeStepSrc();
}


void TurbCalcInvFlux()
{
    UTurbInvFlux * uTurbInvFlux = new UTurbInvFlux();
    uTurbInvFlux->CalcFlux();
    delete uTurbInvFlux;
}

void TurbCalcVisFlux()
{
    UTurbVisFlux * uTurbVisFlux = new UTurbVisFlux();
    uTurbVisFlux->CalcVisFlux();
    delete uTurbVisFlux;
}

void TurbCalcSrcFlux()
{
    UTurbSrcFlux * uTurbSrcFlux = new UTurbSrcFlux();
    uTurbSrcFlux->CalcSrcFlux();
    delete uTurbSrcFlux;
}

void TurbCalcSpectrum()
{
    UTurbSpectrum * uTurbSpectrum = new UTurbSpectrum();
    uTurbSpectrum->CalcSpectrum();
    delete uTurbSpectrum;
}

void TurbCalcDualTimeStepSrc()
{
    //dual time step source
    if ( ctrl.idualtime == 1 )
    {
        UTurbUnsteady * uTurbUnsteady = new UTurbUnsteady();
        uTurbUnsteady->CalcDualTimeSrc();
        delete uTurbUnsteady;
    }
}

void CalcTurbulentViscosity()
{
    UTurbSrcFlux * uTurbSrcFlux = new UTurbSrcFlux();
    uTurbSrcFlux->CalcVist();
    delete uTurbSrcFlux;

}

EndNameSpace