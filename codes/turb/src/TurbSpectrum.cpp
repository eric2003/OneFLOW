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

#include "TurbSpectrum.h"
#include "Iteration.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "TurbLusgs.h"
#include "HXMath.h"
#include "NsIdx.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

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

    //calculate velocity at the cell interface
    Real vnl  = gcom.xfn * ul + gcom.yfn * vl + gcom.zfn * wl - gcom.vfn;
    Real vnr  = gcom.xfn * ur + gcom.yfn * vr + gcom.zfn * wr - gcom.vfn;

    Real absVnl = ABS( vnl );
    Real absVnr = ABS( vnr );

    Real vnRelative = half * ( vnl + vnr );
    Real absVn = ABS( vnRelative ); 

    turbsp.radius1[ ISA ] += gcom.farea * half * absVn;
    turbsp.radius2[ ISA ] += gcom.farea * half * absVn;

    turbsp.matrix1[ ISA ] += gcom.farea * half * ( - vnl - absVn ); //求积分网格在面的右边，网格的法向恰与面的法向相反
    turbsp.matrix2[ ISA ] += gcom.farea * half * (   vnr - absVn );

    //calculate jacobians of the diffusion operator at the cell face i + 1/2
    Real orl     = 1.0 / ( rl + SMALL );
    Real orr     = 1.0 / ( rr + SMALL );
    Real nuLeft  = turbsp.q1[ ISA ];
    Real nuRight = turbsp.q2[ ISA ];

    Real nuet1 = turbcom.visl1 * orl + nuLeft;
    Real nuet2 = turbcom.visl2 * orr + nuRight;

    Real nueff = half * ( nuet1 + nuet2 );
    nueff = ABS( nueff );

    gcom.cvol  = half * ( gcom.cvol1 + gcom.cvol2 );
    Real ovol = one / gcom.cvol;

    Real areaS2 = SQR( gcom.farea );

    turbsp.radius1[ ISA ] +=   turbcom.oreynolds * ovol * ( turbcom.osigma * nueff * areaS2 );
    turbsp.radius2[ ISA ] +=   turbcom.oreynolds * ovol * ( turbcom.osigma * nueff * areaS2 );
            
    turbsp.matrix1[ ISA ] += - turbcom.oreynolds * ovol * ( turbcom.osigma * nueff * areaS2 );
    turbsp.matrix2[ ISA ] += - turbcom.oreynolds * ovol * ( turbcom.osigma * nueff * areaS2 );
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

    //calculate velocity at the cell interface
    Real vnl  = gcom.xfn * ul + gcom.yfn * vl + gcom.zfn * wl - gcom.vfn;
    Real vnr  = gcom.xfn * ur + gcom.yfn * vr + gcom.zfn * wr - gcom.vfn;

    Real absVnl = ABS( vnl );
    Real absVnr = ABS( vnr );

    Real vnRelative = half * ( vnl + vnr );
    Real absVn = ABS( vnRelative ); 

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.radius1[ iEqu ] += gcom.farea * half * absVn;
        turbsp.radius2[ iEqu ] += gcom.farea * half * absVn;

        //即使absVn或者::ABS( vnl )这样的选择也会造成差别，虽然最终对结果的影响还难以预料
        turbsp.matrix1[ iEqu ] += gcom.farea * half * ( - vnl - absVn ); //求积分网格在面的右边，网格的法向恰与面的法向相反
        turbsp.matrix2[ iEqu ] += gcom.farea * half * (   vnr - absVn );
    }

    turbcom.rho  = half * ( rl + rr );
    turbcom.visl = half * ( turbcom.visl1 + turbcom.visl2 );
    turbcom.vist = half * ( turbcom.vist1 + turbcom.vist2 );

    turbcom.CalcSigkw();

    gcom.cvol = half * ( gcom.cvol1 + gcom.cvol2 );
    Real ovol = one / gcom.cvol;

    Real s2 = SQR( gcom.farea );

    turbsp.work[ IKE ] = ( turbcom.visl + turbcom.vist * turbcom.sigk ) / turbcom.rho * s2 * ovol;
    turbsp.work[ IKW ] = ( turbcom.visl + turbcom.vist * turbcom.sigw ) / turbcom.rho * s2 * ovol;

    if ( turbcom.transition_model == ITReGama ) 
    {
        turbsp.work[ ITGama ] = ( turbcom.visl + turbcom.vist / turbcom.trans_df  ) / turbcom.rho * s2 * ovol;
        turbsp.work[ ITRect ] = ( turbcom.visl + turbcom.vist ) * turbcom.trans_dct / turbcom.rho * s2 * ovol;
    }
    /////////////////////////////////////////////////////////////////////////////////////////

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        Real term = turbcom.oreynolds * turbsp.work[ iEqu ];
        turbsp.radius1[ iEqu ] += term;
        turbsp.radius2[ iEqu ] += term;

        turbsp.matrix1[ iEqu ] += - term;
        turbsp.matrix2[ iEqu ] += - term;
    }
}

EndNameSpace