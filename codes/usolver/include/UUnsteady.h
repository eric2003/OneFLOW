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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#pragma once
#include "Unsteady.h"
#include "UnsteadyConvergence.h"

BeginNameSpace( ONEFLOW )

class UUnsteady : public Unsteady
{
public:
    UUnsteady();
public:
    void SetEquationCount( int equationCount );
    void UpdateDualTimeStepResidual();
    void UpdateDualTimeStepSource();
    void StoreOldResidual();
    void PrepareResidual();
    void CalcDualTimeResidual();
    void CalcDualTimeSrc();
    void CalcCellDualTimeResidual();
    void CalcCellDualTimeSrc();
    void CalcUnsteadyCriterion() override;

public:
    // Equation count belongs to the unsteady algorithm state.
    int nEqu;

    // Temporary primitive-state buffers belong to the unsteady algorithm,
    // not to the solver-specific unsteady data object.
    RealField prim, prim1, prim2;

    // Temporary conservative-state buffers belong to the unsteady algorithm,
    // not to the generic unsteady data interface.
    RealField q, q1, q2;
    UnsteadyConvergence convergence;

protected:
    using USDFunc = void( * )( UUnsteady * unst );
    USDFunc srcFun;
    USDFunc criFun;

private:
    RealField res, res1, res2;
    RealField dualtimeRes;
    RealField dualtimeSrc;
};

EndNameSpace
