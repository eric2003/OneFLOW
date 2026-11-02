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

#include "TurbSpectrum.h"
#include "TurbSpectrumPolicy.h"
#include "NsInvFluxPolicy.h"
#include "Iteration.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "TurbLusgs.h"
#include "HXMath.h"
#include "NsIdx.h"
#include "DataBase.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

TurbSpecData turbsp;

TurbSpecData::TurbSpecData()
{
    ;
}

TurbSpecData::~TurbSpecData()
{
    ;
}

void TurbSpecData::Init()
{
    this->nEqu = turbcom.nEqu;
    radius1.resize( nEqu );
    radius2.resize( nEqu );
    matrix1.resize( nEqu );
    matrix2.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );
    work.resize( nEqu );

    int nEqu_ns  = GetDataValue< int >( "nEqu" );

    this->q1_ns.resize( nEqu_ns );
    this->q2_ns.resize( nEqu_ns );
}

TurbSpectrum::TurbSpectrum()
{
    ;
}

TurbSpectrum::~TurbSpectrum()
{
    ;
}

void TurbSpectrum::CalcFaceSpectrum1Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.radius1[ iEqu ] = 0.0;
        turbsp.radius2[ iEqu ] = 0.0;

        turbsp.matrix1[ iEqu ] = 0.0;
        turbsp.matrix2[ iEqu ] = 0.0;
    }

    Real rl = turbsp.q1_ns[ IDX::IR ];
    Real ul = turbsp.q1_ns[ IDX::IU ];
    Real vl = turbsp.q1_ns[ IDX::IV ];
    Real wl = turbsp.q1_ns[ IDX::IW ];

    Real rr = turbsp.q2_ns[ IDX::IR ];
    Real ur = turbsp.q2_ns[ IDX::IU ];
    Real vr = turbsp.q2_ns[ IDX::IV ];
    Real wr = turbsp.q2_ns[ IDX::IW ];

    const Real vnl = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ul, vl, wl );
    const Real vnr = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ur, vr, wr );

    const Real vnRelative = half * ( vnl + vnr );
    const Real absVn = ABS( vnRelative );

    const Real invRadius = TurbInviscidFaceSpectralRadius( gcom.farea, absVn );
    turbsp.radius1[ ISA ] += invRadius;
    turbsp.radius2[ ISA ] += invRadius;

    Real matL, matR;
    TurbInviscidFaceMatrixContrib( gcom.farea, vnl, vnr, absVn, matL, matR );
    turbsp.matrix1[ ISA ] += matL;
    turbsp.matrix2[ ISA ] += matR;

    const Real nuLeft  = turbsp.q1[ ISA ];
    const Real nuRight = turbsp.q2[ ISA ];
    const Real nuet1 = SaEffectiveKinematicViscosity( turbcom.visl1, rl, nuLeft );
    const Real nuet2 = SaEffectiveKinematicViscosity( turbcom.visl2, rr, nuRight );
    Real nueff = ABS( half * ( nuet1 + nuet2 ) );

    gcom.cvol  = half * ( gcom.cvol1 + gcom.cvol2 );
    const Real ovol = one / gcom.cvol;
    const Real areaS2 = SQR( gcom.farea );
    const Real diffRadius = TurbDiffusiveFaceSpectralRadius(
        turbcom.oreynolds, turbcom.osigma, nueff, areaS2, ovol );

    turbsp.radius1[ ISA ] += diffRadius;
    turbsp.radius2[ ISA ] += diffRadius;
    turbsp.matrix1[ ISA ] += - diffRadius;
    turbsp.matrix2[ ISA ] += - diffRadius;
}

void TurbSpectrum::CalcFaceSpectrum2Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.radius1[ iEqu ] = 0.0;
        turbsp.radius2[ iEqu ] = 0.0;

        turbsp.matrix1[ iEqu ] = 0.0;
        turbsp.matrix2[ iEqu ] = 0.0;
    }

    Real rl = turbsp.q1_ns[ IDX::IR ];
    Real ul = turbsp.q1_ns[ IDX::IU ];
    Real vl = turbsp.q1_ns[ IDX::IV ];
    Real wl = turbsp.q1_ns[ IDX::IW ];

    Real rr = turbsp.q2_ns[ IDX::IR ];
    Real ur = turbsp.q2_ns[ IDX::IU ];
    Real vr = turbsp.q2_ns[ IDX::IV ];
    Real wr = turbsp.q2_ns[ IDX::IW ];

    const Real vnl = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ul, vl, wl );
    const Real vnr = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ur, vr, wr );

    const Real vnRelative = half * ( vnl + vnr );
    const Real absVn = ABS( vnRelative );
    const Real invRadius = TurbInviscidFaceSpectralRadius( gcom.farea, absVn );
    Real matL, matR;
    TurbInviscidFaceMatrixContrib( gcom.farea, vnl, vnr, absVn, matL, matR );

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.radius1[ iEqu ] += invRadius;
        turbsp.radius2[ iEqu ] += invRadius;
        turbsp.matrix1[ iEqu ] += matL;
        turbsp.matrix2[ iEqu ] += matR;
    }

    turbcom.rho  = half * ( rl + rr );
    turbcom.visl = half * ( turbcom.visl1 + turbcom.visl2 );
    turbcom.vist = half * ( turbcom.vist1 + turbcom.vist2 );

    turbcom.CalcSigkw();

    gcom.cvol = half * ( gcom.cvol1 + gcom.cvol2 );
    Real ovol = one / gcom.cvol;

    Real s2 = SQR( gcom.farea );

    // work holds nu_eff * A^2 / V before 1/Re scaling
    turbsp.work[ IKE ] = ( turbcom.visl + turbcom.vist * turbcom.sigk ) / turbcom.rho * s2 * ovol;
    turbsp.work[ IKW ] = ( turbcom.visl + turbcom.vist * turbcom.sigw ) / turbcom.rho * s2 * ovol;

    if ( turbcom.transition_model == ITReGama )
    {
        turbsp.work[ ITGama ] = ( turbcom.visl + turbcom.vist / turbcom.trans_df  ) / turbcom.rho * s2 * ovol;
        turbsp.work[ ITRect ] = ( turbcom.visl + turbcom.vist ) * turbcom.trans_dct / turbcom.rho * s2 * ovol;
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        const Real term = TurbDiffusiveFaceSpectralRadius(
            turbcom.oreynolds, one, turbsp.work[ iEqu ], one, one );
        // work already includes A^2/V and nu_eff/rho; pass as nu_eff with unit area/vol
        turbsp.radius1[ iEqu ] += term;
        turbsp.radius2[ iEqu ] += term;
        turbsp.matrix1[ iEqu ] += - term;
        turbsp.matrix2[ iEqu ] += - term;
    }
}

EndNameSpace
