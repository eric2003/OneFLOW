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
#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class Transition
{
public:
    Transition();
    ~Transition();
public:
    Real ce1, ca1, ce2, ca2, dct, df, cct, s1;
public:
    Real ReynoldsNumberBasedOnStrainRate( Real density, Real distanceToWall, Real viscosity,
        Real strainrate, Real reynoldsNumberInflow );
    Real ReynoldsNumberBasedOnDissipation( Real density, Real distanceToWall, Real viscosity,
        Real dissipationRate, Real reynoldsNumberInflow );
    Real TimeScaleInSourceTerm( Real density, Real velocity, Real viscosity, Real reynoldsNumberInflow );
    Real MomentumThickness( Real density, Real velocity, Real viscosity, Real onsetReynoldsOnMomentumThickness,
        Real reynoldsNumberInflow );
    Real FlengthGivenByLangtry( Real localTransitiononsetReynoldsOnMomentumThickness );
    Real HighReynoldsCorrectionOfFlength( Real ReynoldsNumberBasedOnDissipation, Real Flengthori );
    Real ViscosityRatio( Real density, Real viscosity, Real kineticEnergy, Real dissipationRate, Real reynoldsInflow );
    Real TransitionOnsetFunction( Real Rev, Real Rectac, Real RT );
    Real TransitionOnsetMomentumThicknessReynolds( Real Rectabar );
    Real SeparationCorrectionOfIntermittency( Real gmori, Real Rev, Real Rectac, Real RT, Real Fctat );
    Real ControlFunctionFturb( Real RT );
    Real CorrectionOfBlendingFunctionInSST( Real density, Real distance, Real viscosity, Real kineticEnergy, Real reynoldsInflow, Real F1 );
    Real CorrectionOfDestructionInKEquation( Real gmeff, Real Dk );
    Real CorrectionOfProductionInKEquation( Real gmeff, Real Pk );
    Real CalcIntensity( Real velocity, Real kineticEnergy );
    Real EmpiricalCorrelationOfFlamdacta( Real Tu, Real lamdacta );
    Real EmpiricalCorrelationOfRectat( Real Tu, Real Flamdacta );
    Real AccelerationAlongStreamline( Real u, Real v, Real w,
        Real dudx, Real dudy, Real dudz,
        Real dvdx, Real dvdy, Real dvdz,
        Real dwdx, Real dwdy, Real dwdz );
    Real BlendingFunctionOfFctat( Real gm, Real Rew, Real Rectabar, Real vorticity, Real viscosity, Real density, Real absu, Real distance, Real reynoldsInflow );
    Real PressureGradientFunction( Real density, Real momentumThickness, Real viscosity, Real dUds, Real reynoldsInflow );
};

extern Transition turb_trans;


EndNameSpace