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
#include "TurbTrans.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "Ctrl.h"
#include "DataBase.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )
Transition turb_trans;
Transition::Transition()
{
    ;
}

Transition::~Transition()
{
    ;
}

Real Transition::ReynoldsNumberBasedOnStrainRate( Real density, Real distanceToWall, Real viscosity,
                                        Real strainrate, Real reynoldsNumberInflow )
{
    return density * distanceToWall * distanceToWall * strainrate * reynoldsNumberInflow / viscosity;
}

Real Transition::ReynoldsNumberBasedOnDissipation( Real density, Real distanceToWall, Real viscosity,
                                         Real dissipationRate, Real reynoldsNumberInflow )
{
    return density * distanceToWall * distanceToWall * dissipationRate * reynoldsNumberInflow / viscosity;
}

Real Transition::TimeScaleInSourceTerm( Real density, Real velocity, Real viscosity, Real reynoldsNumberInflow )
{
    return 500.0 * viscosity / ( density * velocity * velocity * reynoldsNumberInflow );
}

Real Transition::MomentumThickness( Real density, Real velocity, Real viscosity, Real onsetReynoldsOnMomentumThickness,
                          Real reynoldsNumberInflow )
{
    return onsetReynoldsOnMomentumThickness * viscosity / ( density * velocity * reynoldsNumberInflow );
}

Real Transition::FlengthGivenByLangtry( Real localTransitiononsetReynoldsOnMomentumThickness )
{
    Real Rectabar, Flength;

    Rectabar = localTransitiononsetReynoldsOnMomentumThickness;

    if ( Rectabar < 400.0 ) 
    {
        Flength = 3.98189e+1 - Rectabar * ( 1.1927e-2 + 1.32567e-4 * Rectabar );
    }
    else if ( Rectabar >= 400.0 && Rectabar < 596.0 ) 
    {
       Flength = 2.63404e2 - Rectabar * ( 1.23939 - Rectabar * ( 1.94548e-3 - 1.01695e-6 * Rectabar ) );
    }
    else if ( Rectabar >= 596.0 && Rectabar < 1200. ) 
    {
        Flength = 0.5 - 3.0e-4 * ( Rectabar - 596.0 );
    }
    else
    {
        Flength = 0.3188;
    }
    return Flength;
}

Real Transition::HighReynoldsCorrectionOfFlength( Real ReynoldsNumberBasedOnDissipation, Real Flengthori )
{
    Real Rw = ReynoldsNumberBasedOnDissipation / 500.0;
    Real Fsublayer = exp( -( 6.25 * Rw * Rw ) );
    return Flengthori - ( Flengthori - 40.0 ) * Fsublayer;
}

Real Transition::ViscosityRatio( Real density, Real viscosity, Real kineticEnergy, Real dissipationRate, Real reynoldsInflow )
{
    return density * kineticEnergy * reynoldsInflow / ( viscosity * dissipationRate );
}

Real Transition::TransitionOnsetFunction( Real Rev, Real Rectac, Real RT )
{
    Real Fonset1 = Rev / ( 2.193 * Rectac );
    Real Fonset2 = MIN( MAX( Fonset1, pow( Fonset1, 4.0 ) ), 2.0 );
    Real Fonset3 = MAX( 1.0 - RT * RT * RT / 15.625 , 0.0 );
       
    return MAX( Fonset2 - Fonset3, 0.0 );
}

Real Transition::TransitionOnsetMomentumThicknessReynolds( Real Rectabar )
{
    if ( Rectabar > 1870.0 ) 
    {
        return Rectabar - 5.9311e2 - 0.482 * ( Rectabar - 1870.0 );
    }
    else
    {
        return -3.96035 + Rectabar * ( 1.0120656 - Rectabar * ( 8.6823e-4 - Rectabar * ( 6.96506e-7 - 1.74105e-10 * Rectabar ) ) );
    }
};

Real Transition::SeparationCorrectionOfIntermittency( Real gmori, Real Rev, Real Rectac, Real RT, Real Fctat )
{
    Real Freattach = exp( - pow( RT / 20.0, 4.0 ) );

    Real gmsep = Fctat * MIN( 2.0, s1 * Freattach * MAX( 0.0, Rev / ( 3.235 * Rectac ) - 1.0 ) );
    
    return MAX( gmori, gmsep );
}

Real Transition::CalcIntensity( Real velocity, Real kineticEnergy )
{
    return MAX( 0.027, 100.0 * sqrt( 2.0 * kineticEnergy / 3.0 ) / velocity );
}

Real Transition::ControlFunctionFturb( Real RT )
{
    return exp( -pow( RT / 4.0, 4.0 ) );
}

Real Transition::CorrectionOfBlendingFunctionInSST( Real density, Real distance, Real viscosity, Real kineticEnergy, Real reynoldsInflow, Real F1 )
{
    Real ry = reynoldsInflow * density * distance * sqrt( kineticEnergy ) / viscosity;
    
    return MAX( F1, exp( -pow( ry / 120.0, 8.0 ) ) );
}

Real Transition::CorrectionOfDestructionInKEquation( Real gmeff, Real Dk )
{
    return Dk * MIN( 1.0, MAX( 0.1, gmeff ) );
}

Real Transition::CorrectionOfProductionInKEquation( Real gmeff, Real Pk )
{
    return gmeff * Pk;
}

Real Transition::EmpiricalCorrelationOfFlamdacta( Real Tu, Real lamdacta )
{
    if ( lamdacta > 0.0 ) 
    {
        return 1.0 + 0.275 * ( 1.0 - exp(-35.0 * lamdacta ) ) *  exp(-2.0 * Tu );
    }
    else if ( lamdacta < 0.0 ) 
    {
        return 1.0 + lamdacta * ( 12.986 + lamdacta * ( 123.66 + 405.689 * lamdacta ) ) * exp(-pow( Tu / 1.5, 1.5 ) );
    }
    else
    {
        return 1.0;
    }
}

Real Transition::EmpiricalCorrelationOfRectat( Real Tu, Real Flamdacta )
{
    Real Rectat;
    if ( Tu > 1.3 ) 
    {
        Rectat = 331.5 * Flamdacta * pow( Tu - 0.5658,-0.671 );
    }
    else
    {
        Rectat = Flamdacta * ( 1173.51 - 589.428 * Tu + 0.2196 / ( Tu * Tu ) );
    }
       
    return MAX( 20.0, Rectat );
}

Real Transition::AccelerationAlongStreamline( Real u   , Real v   , Real w   , 
                                    Real dudx, Real dudy, Real dudz,
                                    Real dvdx, Real dvdy, Real dvdz,
                                    Real dwdx, Real dwdy, Real dwdz )
{
    Real absu = MAX( 1.0e-20, DIST( u, v, w ) );
    
    Real uu = u / absu;
    Real vv = v / absu;
    Real ww = w / absu;
    
    Real dUx = uu * dudx + vv * dvdx + ww * dwdx;
    Real dUy = uu * dudy + vv * dvdy + ww * dwdy;
    Real dUz = uu * dudz + vv * dvdz + ww * dwdz;
    
    return uu * dUx + vv * dUy + ww * dUz;
}

Real Transition::BlendingFunctionOfFctat( Real gm, Real Rew, Real Rectabar, Real vorticity, Real viscosity, Real density, Real absu, Real distance, Real reynoldsInflow )
{
    Real ctaBL = Rectabar * viscosity / MAX( reynoldsInflow * density * absu, 1.0e-20 );
    
    Real deltaBL = 7.5 * ctaBL;
    
    Real delta = 50.0 * vorticity * distance * deltaBL / MAX( absu, 1.0e-20 );
    
    Real Fwake = exp( - ( Rew / 1.0e5 ) * ( Rew / 1.0e5 ) );

    Real tmp1 = Fwake * exp( - pow( distance / delta, 4.0 ) );
    Real tmp2 = 1.0 - ( ( ce2 * gm - 1.0 ) / ( ce2 - 1.0 ) ) * ( ( ce2 * gm - 1.0 ) / ( ce2 - 1.0 ) );
       
    return MIN( 1.0, MAX( tmp1, tmp2 ) );
}

Real Transition::PressureGradientFunction( Real density, Real momentumThickness, Real viscosity, Real dUds, Real reynoldsInflow )
{        
    Real lamdacta = density * momentumThickness * momentumThickness * dUds * reynoldsInflow / viscosity;
    
    return MAX( -0.1, MIN( 0.1, lamdacta ) );
}


EndNameSpace