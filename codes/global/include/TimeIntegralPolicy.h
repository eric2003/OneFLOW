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

BeginNameSpace( ONEFLOW )

// @doc-algorithm time-integration
// @doc-stage rk-stage-weight
// @doc-pure
// @doc-formula a_s = 1 / (s + 1),  s = 0, ..., N-1
// OneFLOW multi-stage coefficient schedule (not a classical Butcher tableau).
// See doc/source/numerical/time_integral_rk_stage_weight.rst.
inline Real RungeKuttaStageWeight( int stageIndex )
{
    return 1.0 / static_cast< Real >( stageIndex + 1 );
}

EndNameSpace
