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

#include "TurbRhs.h"
#include "Ctrl.h"
#include "UTurbInvFlux.h"
#include "UTurbVisFlux.h"
#include "UTurbSrcFlux.h"
#include "UTurbSpectrum.h"
#include "UTurbUnsteady.h"
#include "UTurbBcSolver.h"
#include <memory> // Added for std::make_unique

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
    UTurbBcSolver uTurbBcSolver;
    uTurbBcSolver.Init();
    uTurbBcSolver.CalcBc();
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
    UTurbInvFlux uTurbInvFlux;
    uTurbInvFlux.CalcFlux();
}

void TurbCalcVisFlux()
{
    UTurbVisFlux uTurbVisFlux;
    uTurbVisFlux.CalcVisFlux();
}

void TurbCalcSrcFlux()
{
    UTurbSrcFlux uTurbSrcFlux;
    uTurbSrcFlux.CalcSrcFlux();
}

void TurbCalcSpectrum()
{
    UTurbSpectrum uTurbSpectrum;
    uTurbSpectrum.CalcSpectrum();
}

void TurbCalcDualTimeStepSrc()
{
    //dual time step source
    if ( ctrl.idualtime == 1 )
    {
        UTurbUnsteady uTurbUnsteady;
        uTurbUnsteady.CalcDualTimeSrc();
    }
}

void CalcTurbulentViscosity()
{
    UTurbSrcFlux uTurbSrcFlux;
    uTurbSrcFlux.CalcVist();
}

EndNameSpace
