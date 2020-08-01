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

#include "INsLusgs.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ILusgsData inslu;

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
        cout << "Navier Stokes residual reduced by " << inslu.dmax << " with " << inslu.numberOfRealSweeps << " Sweeps\n";
    }
}

void INsLusgs::ZeroFluxIncrement()
{
    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.rhs0[ iEqu ] = 0.0;
    }
}

void INsLusgs::AddViscousTerm()
{
    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.rhs0[ iEqu ] -= inslu.visrad * inslu.dqj[ iEqu ];
    }
}

void INsLusgs::AddFluxIncrement()
{
    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.rhs0[ iEqu ] += inslu.dfj[ iEqu ];
    }
}

void INsLusgs::AddFluxIncrement( const Real & coef )
{
    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.rhs0[ iEqu ] += coef * inslu.dfj[ iEqu ];
    }
}

void INsLusgs::GetFluxIncrement( int signOfMatrix )
{
    this->GetStandardFluxIncrement( signOfMatrix );
}

void INsLusgs::CalcFaceEigenValue( RealField & prim )
{
    //这里输入的应该是作用单元界面上的值
    Real & rm  = prim[ IIDX::IIR ];
    Real & um  = prim[ IIDX::IIU ];
    Real & vm  = prim[ IIDX::IIV ];
    Real & wm  = prim[ IIDX::IIW ];
    Real & pm  = prim[ IIDX::IIP ];

    Real c2 = ABS( inslu.gama * pm / rm );
    Real cm = sqrt( c2 );

    Real vn_fluid = gcom.xfn * um + gcom.yfn * vm + gcom.zfn * wm;
    Real vn_rel   = vn_fluid - gcom.vfn;

    //Real lmd1 = vn_rel;
    //Real lmd2 = vn_rel + cm;
    //Real lmd3 = vn_rel - cm;

    Real max_eigen = ABS( vn_rel ) + cm;

    inslu.lmdOnFace1 = max_eigen;
    inslu.lmdOnFace2 = max_eigen;
    inslu.lmdOnFace3 = max_eigen;
}

void INsLusgs::GetStandardFluxIncrement( int signOfMatrix )
{
    this->CalcFaceEigenValue( inslu.primF );

    Real & rm  = inslu.primj[ IIDX::IIR ];
    Real & um  = inslu.primj[ IIDX::IIU ];
    Real & vm  = inslu.primj[ IIDX::IIV ];
    Real & wm  = inslu.primj[ IIDX::IIW ];
    Real & pm  = inslu.primj[ IIDX::IIP ];

    Real c2 = ABS( inslu.gama * pm / rm );
    Real cm = sqrt( c2 );

    Real vn_fluid = gcom.xfn * um + gcom.yfn * vm + gcom.zfn * wm;
    Real vn_rel   = vn_fluid - gcom.vfn;

    Real lmd1 = vn_rel;
    Real lmd2 = vn_rel + cm;
    Real lmd3 = vn_rel - cm;

    lmd1 = half * ( lmd1 + signOfMatrix * inslu.lmdOnFace1 );
    lmd2 = half * ( lmd2 + signOfMatrix * inslu.lmdOnFace2 );
    lmd3 = half * ( lmd3 + signOfMatrix * inslu.lmdOnFace3 );

    Real x1 = ( lmd1 + lmd1 - lmd2 - lmd3 ) / ( c2 + c2 );
    Real x2 = ( lmd2 - lmd3 ) / ( cm + cm );

    Real dc =   vn_fluid * inslu.dqj[ IIDX::IIR  ]
              - gcom.xfn * inslu.dqj[ IIDX::IIRU ]
              - gcom.yfn * inslu.dqj[ IIDX::IIRV ]
              - gcom.zfn * inslu.dqj[ IIDX::IIRW ];
    Real c2dc = c2 * dc;

    Real dh, hm;

    ONEFLOW::CalcIDH( inslu.primj, inslu.gama, inslu.dqj, dh, hm );

    Real term1 =  dh   * x1 + dc * x2;
    Real term2 =  c2dc * x1 + dh * x2;

    inslu.dfj[ IIDX::IIR  ] = lmd1 * inslu.dqj[ IIDX::IIR  ] -      term1                    ;
    inslu.dfj[ IIDX::IIRU ] = lmd1 * inslu.dqj[ IIDX::IIRU ] - um * term1 + gcom.xfn * term2;
    inslu.dfj[ IIDX::IIRV ] = lmd1 * inslu.dqj[ IIDX::IIRV ] - vm * term1 + gcom.yfn * term2;
    inslu.dfj[ IIDX::IIRW ] = lmd1 * inslu.dqj[ IIDX::IIRW ] - wm * term1 + gcom.zfn * term2;
   // inslu.dfj[ IIDX::IIRE ] = lmd1 * inslu.dqj[ IIDX::IIRE ] - hm * term1 + vn_fluid * term2;

    for ( int iEqu = inslu.nBEqu; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.dfj[ iEqu ] = lmd1 * inslu.dqj[ iEqu ] - inslu.primj[ iEqu ] * term1;
    }

    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.dfj[ iEqu ] *= gcom.farea;
    }
}


void INsLusgs::InitializeSweep( int iSweep )
{
    inslu.norm = 0.0;
}

bool INsLusgs::UpdateSweep( int iSweep )
{
    inslu.numberOfRealSweeps = iSweep + 1;
    if ( iSweep == 0 )
    {
        inslu.norm0 = inslu.norm;
        inslu.dmax  = 1.0;
    }
    else
    {
        inslu.dmax = sqrt( inslu.norm / ( inslu.norm0 + SMALL ) );
    }

    if ( inslu.dmax < inslu.tol )
    {
        return true;
    }
    return false;
}

void INsLusgs::CalcLowerChange()
{
    if ( inslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] = inslu.dqi[ iEqu ] - inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] /=  inslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqi[ iEqu ] = inslu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.drhs[ iEqu ] += inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqSweep =  inslu.dqi[ iEqu ] - inslu.dqi0[ iEqu ];
            inslu.norm   += SQR( inslu.dqSweep );
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] = inslu.rhs[ iEqu ] - inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] /=  inslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqi[ iEqu ] = inslu.tmp[ iEqu ];
        }
    }
}

void INsLusgs::CalcUpperChange()
{
    if ( inslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] = inslu.dqi[ iEqu ] - inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] /=  inslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqi[ iEqu ] = inslu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {            
            inslu.drhs[ iEqu ] += inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqSweep     = inslu.dqi[ iEqu ] - inslu.dqi0[ iEqu ];
            inslu.norm       += SQR( inslu.dqSweep );
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] = - inslu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.tmp[ iEqu ] /=  inslu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
        {
            inslu.dqi[ iEqu ] += inslu.tmp[ iEqu ];
        }
    }
}


bool INsLusgs::IsOversetCell()
{
    return ( gcom.blank <= 0 );
}

void INsLusgs::ZeroOversetCell()
{
    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.dqi[ iEqu ] = 0.0;
    }

    for ( int iEqu = 0; iEqu < inslu.nEqu; ++ iEqu )
    {
        inslu.drhs[ iEqu ] = 0.0;
    }
}

void CalcIDH(RealField & prim, Real & gama, RealField & dq, Real & dh, Real & totalEnthalpy)
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