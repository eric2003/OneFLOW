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

#include "NsLusgs.h"
#include "NsLusgsPolicy.h"
#include "NsInvFluxPolicy.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

LusgsData nslu;

void LusgsData::Init()
{
    nEqu = nscom.nEqu;
    nBEqu = nEqu;

    radius.resize( nEqu );
    dqj.resize( nEqu );
    dqi.resize( nEqu );
    dqi0.resize( nEqu );
    primj.resize( nEqu );
    primF.resize( nEqu );
    rhs0.resize( nEqu );
    dfj.resize( nEqu );
    drhs.resize( nEqu );
    rhs.resize( nEqu );
    tmp.resize( nEqu );
}

void NsLusgs::InitializeSub()
{
}

void NsLusgs::DumpSweepInformation()
{
    int pid = ONEFLOW::Parallel::GetPid();
    if ( pid == ONEFLOW::Parallel::GetServerid() )
    {
        std::cout << "Navier Stokes residual reduced by " << nslu.dmax << " with " << nslu.numberOfRealSweeps << " Sweeps\n";
    }
}

void NsLusgs::ZeroFluxIncrement()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] = 0.0;
    }
}

void NsLusgs::AddViscousTerm()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] -= nslu.visrad * nslu.dqj[ iEqu ];
    }
}

void NsLusgs::AddFluxIncrement()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] += nslu.dfj[ iEqu ];
    }
}

void NsLusgs::AddFluxIncrement( const Real & coef )
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] += coef * nslu.dfj[ iEqu ];
    }
}

void NsLusgs::GetFluxIncrement( int signOfMatrix )
{
    this->GetStandardFluxIncrement( signOfMatrix );
}

void NsLusgs::CalcFaceEigenValue( RealField & prim )
{
    // Face state on the LUSGS interface unit.
    Real & rm  = prim[ IDX::IR ];
    Real & um  = prim[ IDX::IU ];
    Real & vm  = prim[ IDX::IV ];
    Real & wm  = prim[ IDX::IW ];
    Real & pm  = prim[ IDX::IP ];

    const Real cm = IdealGasSoundSpeedFromPressure( rm, pm, nslu.gama );
    const Real vn_rel = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, um, vm, wm );

    const Real max_eigen = ABS( vn_rel ) + cm;

    nslu.lmdOnFace1 = max_eigen;
    nslu.lmdOnFace2 = max_eigen;
    nslu.lmdOnFace3 = max_eigen;
}

void NsLusgs::GetStandardFluxIncrement( int signOfMatrix )
{
    this->CalcFaceEigenValue( nslu.primF );

    Real & rm  = nslu.primj[ IDX::IR ];
    Real & um  = nslu.primj[ IDX::IU ];
    Real & vm  = nslu.primj[ IDX::IV ];
    Real & wm  = nslu.primj[ IDX::IW ];
    Real & pm  = nslu.primj[ IDX::IP ];

    const Real cm = IdealGasSoundSpeedFromPressure( rm, pm, nslu.gama );
    const Real c2 = cm * cm;
    const Real vn_fluid = gcom.xfn * um + gcom.yfn * vm + gcom.zfn * wm;
    const Real vn_rel = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, um, vm, wm );

    Real lmd1, lmd2, lmd3;
    LusgsSplitFaceEigenvalues(
        vn_rel, cm, static_cast< Real >( signOfMatrix ),
        nslu.lmdOnFace1, nslu.lmdOnFace2, nslu.lmdOnFace3,
        lmd1, lmd2, lmd3 );

    Real x1, x2;
    LusgsLinearizationWeights( lmd1, lmd2, lmd3, c2, cm, x1, x2 );

    const Real dc = RoeContinuityNormalJump(
        vn_fluid, gcom.xfn, gcom.yfn, gcom.zfn,
        nslu.dqj[ IDX::IR ], nslu.dqj[ IDX::IRU ],
        nslu.dqj[ IDX::IRV ], nslu.dqj[ IDX::IRW ] );

    Real dh, hm;
    ONEFLOW::CalcDH( nslu.primj, nslu.gama, nslu.dqj, dh, hm );

    Real term1, term2;
    RoeDissipationCouplingTerms( dh, dc, c2, x1, x2, term1, term2 );

    LusgsFluxIncrementFiveEq(
        lmd1, term1, term2,
        um, vm, wm, hm, vn_fluid,
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.farea,
        nslu.dqj[ IDX::IR ], nslu.dqj[ IDX::IRU ], nslu.dqj[ IDX::IRV ],
        nslu.dqj[ IDX::IRW ], nslu.dqj[ IDX::IRE ],
        nslu.dfj[ IDX::IR ], nslu.dfj[ IDX::IRU ], nslu.dfj[ IDX::IRV ],
        nslu.dfj[ IDX::IRW ], nslu.dfj[ IDX::IRE ] );

    for ( int iEqu = nslu.nBEqu; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dfj[ iEqu ] = ( lmd1 * nslu.dqj[ iEqu ] - nslu.primj[ iEqu ] * term1 )
            * gcom.farea;
    }
}


void NsLusgs::InitializeSweep( int iSweep )
{
    nslu.norm = 0.0;
}

bool NsLusgs::UpdateSweep( int iSweep )
{
    nslu.numberOfRealSweeps = iSweep + 1;
    if ( iSweep == 0 )
    {
        nslu.norm0 = nslu.norm;
        nslu.dmax  = 1.0;
    }
    else
    {
        nslu.dmax = sqrt( nslu.norm / ( nslu.norm0 + SMALL ) );
    }

    if ( nslu.dmax < nslu.tol )
    {
        return true;
    }
    return false;
}

void NsLusgs::CalcLowerChange()
{
    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] = nslu.dqi[ iEqu ] - nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] /=  nslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi[ iEqu ] = nslu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.drhs[ iEqu ] += nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqSweep =  nslu.dqi[ iEqu ] - nslu.dqi0[ iEqu ];
            nslu.norm   += SQR( nslu.dqSweep );
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] = nslu.rhs[ iEqu ] - nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] /=  nslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi[ iEqu ] = nslu.tmp[ iEqu ];
        }
    }
}

void NsLusgs::CalcUpperChange()
{
    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] = nslu.dqi[ iEqu ] - nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] /=  nslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi[ iEqu ] = nslu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {            
            nslu.drhs[ iEqu ] += nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqSweep     = nslu.dqi[ iEqu ] - nslu.dqi0[ iEqu ];
            nslu.norm       += SQR( nslu.dqSweep );
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] = - nslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.tmp[ iEqu ] /=  nslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi[ iEqu ] += nslu.tmp[ iEqu ];
        }
    }
}


bool NsLusgs::IsOversetCell()
{
    return ( gcom.blank <= 0 );
}

void NsLusgs::ZeroOversetCell()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dqi[ iEqu ] = 0.0;
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.drhs[ iEqu ] = 0.0;
    }
}

void CalcDH( RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy )
{
    Real & density  = prim[ IDX::IR ];
    Real & um       = prim[ IDX::IU ];
    Real & vm       = prim[ IDX::IV ];
    Real & wm       = prim[ IDX::IW ];
    Real & pressure = prim[ IDX::IP ];

    const Real v2 = ONEFLOW::SQR( um, vm, wm );
    dh = RoeTotalEnthalpyJump(
        um, vm, wm, gama,
        dq[ IDX::IR ], dq[ IDX::IRU ], dq[ IDX::IRV ], dq[ IDX::IRW ], dq[ IDX::IRE ] );
    totalEnthalpy = TotalEnthalpy(
        IdealGasStaticEnthalpy( density, pressure, gama ), v2 );
}


EndNameSpace
