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

#include "NsLusgs.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

LusgsData nslu;

LusgsData::LusgsData()
{
    ;
}

LusgsData::~LusgsData()
{
    ;
}

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

NsLusgs::NsLusgs()
{
}

NsLusgs::~NsLusgs()
{
}

void NsLusgs::InitializeSub()
{
}

void NsLusgs::DumpSweepInformation()
{
    int pid = ONEFLOW::Parallel::GetPid();
    if ( pid == ONEFLOW::Parallel::GetServerid() )
    {
        cout << "Navier Stokes residual reduced by " << nslu.dmax << " with " << nslu.numberOfRealSweeps << " Sweeps\n";
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
    //这里输入的应该是作用单元界面上的值
    Real & rm  = prim[ IDX::IR ];
    Real & um  = prim[ IDX::IU ];
    Real & vm  = prim[ IDX::IV ];
    Real & wm  = prim[ IDX::IW ];
    Real & pm  = prim[ IDX::IP ];

    Real c2 = ABS( nslu.gama * pm / rm );
    Real cm = sqrt( c2 );

    Real vn_fluid = gcom.xfn * um + gcom.yfn * vm + gcom.zfn * wm;
    Real vn_rel   = vn_fluid - gcom.vfn;

    //Real lmd1 = vn_rel;
    //Real lmd2 = vn_rel + cm;
    //Real lmd3 = vn_rel - cm;

    Real max_eigen = ABS( vn_rel ) + cm;

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

    Real c2 = ABS( nslu.gama * pm / rm );
    Real cm = sqrt( c2 );

    Real vn_fluid = gcom.xfn * um + gcom.yfn * vm + gcom.zfn * wm;
    Real vn_rel   = vn_fluid - gcom.vfn;

    Real lmd1 = vn_rel;
    Real lmd2 = vn_rel + cm;
    Real lmd3 = vn_rel - cm;

    lmd1 = half * ( lmd1 + signOfMatrix * nslu.lmdOnFace1 );
    lmd2 = half * ( lmd2 + signOfMatrix * nslu.lmdOnFace2 );
    lmd3 = half * ( lmd3 + signOfMatrix * nslu.lmdOnFace3 );

    Real x1 = ( lmd1 + lmd1 - lmd2 - lmd3 ) / ( c2 + c2 );
    Real x2 = ( lmd2 - lmd3 ) / ( cm + cm );

    Real dc =   vn_fluid * nslu.dqj[ IDX::IR  ]
              - gcom.xfn * nslu.dqj[ IDX::IRU ]
              - gcom.yfn * nslu.dqj[ IDX::IRV ]
              - gcom.zfn * nslu.dqj[ IDX::IRW ];
    Real c2dc = c2 * dc;

    Real dh, hm;

    ONEFLOW::CalcDH( nslu.primj, nslu.gama, nslu.dqj, dh, hm );

    Real term1 =  dh   * x1 + dc * x2;
    Real term2 =  c2dc * x1 + dh * x2;

    nslu.dfj[ IDX::IR  ] = lmd1 * nslu.dqj[ IDX::IR  ] -      term1                    ;
    nslu.dfj[ IDX::IRU ] = lmd1 * nslu.dqj[ IDX::IRU ] - um * term1 + gcom.xfn * term2;
    nslu.dfj[ IDX::IRV ] = lmd1 * nslu.dqj[ IDX::IRV ] - vm * term1 + gcom.yfn * term2;
    nslu.dfj[ IDX::IRW ] = lmd1 * nslu.dqj[ IDX::IRW ] - wm * term1 + gcom.zfn * term2;
    nslu.dfj[ IDX::IRE ] = lmd1 * nslu.dqj[ IDX::IRE ] - hm * term1 + vn_fluid * term2;

    for ( int iEqu = nslu.nBEqu; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dfj[ iEqu ] = lmd1 * nslu.dqj[ iEqu ] - nslu.primj[ iEqu ] * term1;
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dfj[ iEqu ] *= gcom.farea;
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
    Real v2, ae, af;
    Real enthalpy;

    Real & density  = prim[ IDX::IR ];
    Real & um       = prim[ IDX::IU ];
    Real & vm       = prim[ IDX::IV ];
    Real & wm       = prim[ IDX::IW ];
    Real & pressure = prim[ IDX::IP ];

    v2 = ONEFLOW::SQR( um, vm, wm );

    ae = gama - one;
    af = half * ae * v2;
    dh = - ae * ( um * dq[ IDX::IRU ] + vm * dq[ IDX::IRV ] + wm * dq[ IDX::IRW ] - dq[ IDX::IRE ] );

    enthalpy = ( gama / ( gama - one ) ) * ( pressure / density );

    dh += af * dq[ IDX::IR  ];
    totalEnthalpy = enthalpy + half * v2;
}


EndNameSpace