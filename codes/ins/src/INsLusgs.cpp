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

#include "INsLusgs.h"
#include "INsCom.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Parallel.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ILusgsData nslu;

ILusgsData::ILusgsData()
{
    ;
}

ILusgsData::~ILusgsData()
{
    ;
}

void ILusgsData::Init()
{
    nEqu = inscom.nEqu;
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

INsLusgs::INsLusgs()
{
}

INsLusgs::~INsLusgs()
{
}

void INsLusgs::InitializeSub()
{
}

void INsLusgs::DumpSweepInformation()
{
    int pid = ONEFLOW::Parallel::GetPid();
    if ( pid == ONEFLOW::Parallel::GetServerid() )
    {
        cout << "Navier Stokes residual reduced by " << nslu.dmax << " with " << nslu.numberOfRealSweeps << " Sweeps\n";
    }
}

void INsLusgs::ZeroFluxIncrement()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] = 0.0;
    }
}

void INsLusgs::AddViscousTerm()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] -= nslu.visrad * nslu.dqj[ iEqu ];
    }
}

void INsLusgs::AddFluxIncrement()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] += nslu.dfj[ iEqu ];
    }
}

void INsLusgs::AddFluxIncrement( const Real & coef )
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] += coef * nslu.dfj[ iEqu ];
    }
}

void INsLusgs::GetFluxIncrement( int signOfMatrix )
{
    this->GetStandardFluxIncrement( signOfMatrix );
}

void INsLusgs::CmpFaceEigenValue( RealField & prim )
{
    //这里输入的应该是作用单元界面上的值
    Real & rm  = prim[ IIDX::IIR ];
    Real & um  = prim[ IIDX::IIU ];
    Real & vm  = prim[ IIDX::IIV ];
    Real & wm  = prim[ IIDX::IIW ];
    Real & pm  = prim[ IIDX::IIP ];

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

void INsLusgs::GetStandardFluxIncrement( int signOfMatrix )
{
    this->CmpFaceEigenValue( nslu.primF );

    Real & rm  = nslu.primj[ IIDX::IIR ];
    Real & um  = nslu.primj[ IIDX::IIU ];
    Real & vm  = nslu.primj[ IIDX::IIV ];
    Real & wm  = nslu.primj[ IIDX::IIW ];
    Real & pm  = nslu.primj[ IIDX::IIP ];

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

    Real dc =   vn_fluid * nslu.dqj[ IIDX::IIR  ]
              - gcom.xfn * nslu.dqj[ IIDX::IIRU ]
              - gcom.yfn * nslu.dqj[ IIDX::IIRV ]
              - gcom.zfn * nslu.dqj[ IIDX::IIRW ];
    Real c2dc = c2 * dc;

    Real dh, hm;

    ONEFLOW::CmpIDH( nslu.primj, nslu.gama, nslu.dqj, dh, hm );

    Real term1 =  dh   * x1 + dc * x2;
    Real term2 =  c2dc * x1 + dh * x2;

    nslu.dfj[ IIDX::IIR  ] = lmd1 * nslu.dqj[ IIDX::IIR  ] -      term1                    ;
    nslu.dfj[ IIDX::IIRU ] = lmd1 * nslu.dqj[ IIDX::IIRU ] - um * term1 + gcom.xfn * term2;
    nslu.dfj[ IIDX::IIRV ] = lmd1 * nslu.dqj[ IIDX::IIRV ] - vm * term1 + gcom.yfn * term2;
    nslu.dfj[ IIDX::IIRW ] = lmd1 * nslu.dqj[ IIDX::IIRW ] - wm * term1 + gcom.zfn * term2;
    nslu.dfj[ IIDX::IIRE ] = lmd1 * nslu.dqj[ IIDX::IIRE ] - hm * term1 + vn_fluid * term2;

    for ( int iEqu = nslu.nBEqu; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dfj[ iEqu ] = lmd1 * nslu.dqj[ iEqu ] - nslu.primj[ iEqu ] * term1;
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dfj[ iEqu ] *= gcom.farea;
    }
}


void INsLusgs::InitializeSweep( int iSweep )
{
    nslu.norm = 0.0;
}

bool INsLusgs::UpdateSweep( int iSweep )
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

void INsLusgs::CmpLowerChange()
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

void INsLusgs::CmpUpperChange()
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


bool INsLusgs::IsOversetCell()
{
    return ( gcom.blank <= 0 );
}

void INsLusgs::ZeroOversetCell()
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

void CmpIDH(RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy)
{
	Real v2, ae, af;
	Real enthalpy;

	Real & density = prim[IIDX::IIR];
	Real & um = prim[IIDX::IIU];
	Real & vm = prim[IIDX::IIV];
	Real & wm = prim[IIDX::IIW];
	Real & pressure = prim[IIDX::IIP];

	v2 = ONEFLOW::SQR(um, vm, wm);

	ae = gama - one;
	af = half * ae * v2;
	//dh = -ae * (um * dq[IIDX::IIRU] + vm * dq[IIDX::IIRV] + wm * dq[IIDX::IIRW] - dq[IIDX::IIRE]);

	enthalpy = (gama / (gama - one)) * (pressure / density);

	dh += af * dq[IIDX::IIR];
	totalEnthalpy = enthalpy + half * v2;
}




EndNameSpace