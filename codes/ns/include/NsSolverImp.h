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
#include "HXClone.h"

BeginNameSpace( ONEFLOW )

DEFINE_DATA_CLASS( NsInitFinal );
DEFINE_DATA_CLASS( NsVisual );
DEFINE_DATA_CLASS( NsCmpTimeStep );
DEFINE_DATA_CLASS( NsUpdateResiduals );
DEFINE_DATA_CLASS( NsImplicitMethod );
DEFINE_DATA_CLASS( NsPostprocess );
DEFINE_DATA_CLASS( NsFinalPostprocess );
DEFINE_DATA_CLASS( NsInitSolver );
DEFINE_DATA_CLASS( NsCmpBoundary );
DEFINE_DATA_CLASS( DumpHeatFluxCoeff );

void RegisterNsFunc();

class SolverRegData;
SolverRegData * GetNsReg();

EndNameSpace