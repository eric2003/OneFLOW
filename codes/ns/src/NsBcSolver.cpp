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

#include "NsBcSolver.h"
#include "BcData.h"
#include "NsCom.h"
#include "UCom.h"
#include "NsCtrl.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Stop.h"
#include "Boundary.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

NsBcSolver::NsBcSolver()
{
    ;
}

NsBcSolver::~NsBcSolver()
{
    ;
}

void NsBcSolver::SetBc()
{
    this->bcPointer = & NsBcSolver::NoBc;

    this->updateFlag = true;

    if ( ug.bctype < 0 )
    {
        this->bcPointer  = & NsBcSolver::InterfaceBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::OVERSET )
    {
        this->bcPointer  = & NsBcSolver::OversetBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::EXTRAPOLATION )
    {
        this->bcPointer = & NsBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::SOLID_SURFACE )
    {
        this->SetSolidSurfaceBc();
    }
    else if ( ug.bctype == BC::SYMMETRY )
    {
        this->bcPointer = & NsBcSolver::SymmetryBc;
    }
    else if ( ug.bctype == BC::FARFIELD )
    {
        this->bcPointer = & NsBcSolver::FarFieldBc;
    }
    else if ( ug.bctype == BC::INFLOW )
    {
        this->bcPointer = & NsBcSolver::InFlowBc;
    }
    else if ( ug.bctype == BC::OUTFLOW )
    {
        this->bcPointer = & NsBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE )
    {
        this->bcPointer = & NsBcSolver::PoleBc;
    }
    else if ( ug.bctype == BC::GENERIC_2 )
    {
        //userDefined
        this->bcPointer = & NsBcSolver::UserDefinedBc;
    }
    else if ( ug.bctype == BC::PERIODIC ) 
    {
        this->bcPointer = & NsBcSolver::PeriodicBc;
    }
    else
    {
        cout << "Error : Illegal BCtype ID " << ug.bctype << endl;
        Stop("");
    }
}

void NsBcSolver::SetSolidSurfaceBc()
{
    if ( vis_model.vismodel == INVISCID )
    {
        this->bcPointer = & NsBcSolver::SymmetryBc;
    }
    else if ( ns_ctrl.isowallbc == 0 )
    {
        //viscous adiabatic wall
        this->bcPointer = & NsBcSolver::AdiabaticVisWallBc;
    }
    else
    {
        //viscous iso-thermal wall
        this->bcPointer = & NsBcSolver::IsothermalVisWallBc;
    }
}

void NsBcSolver::CalcFaceBc()
{
    ( this->* bcPointer )();
}

void NsBcSolver::InFlowBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.inflow[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.inflow[ iEqu ];
    }
}

void NsBcSolver::OutFlowBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
}

void NsBcSolver::PoleBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
}

void NsBcSolver::FarFieldBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    //inner point
    Real rin, uin, vin, win, pin;
    ONEFLOW::Extract( nscom.prims1, rin, uin, vin, win, pin );

    gcom.xfn *= nscom.faceOuterNormal;
    gcom.yfn *= nscom.faceOuterNormal;
    gcom.zfn *= nscom.faceOuterNormal;

    Real rref = nscom.inflow[ IDX::IR ];
    Real uref = nscom.inflow[ IDX::IU ];
    Real vref = nscom.inflow[ IDX::IV ];
    Real wref = nscom.inflow[ IDX::IW ];
    Real pref = nscom.inflow[ IDX::IP ];

    Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
    Real vnin  = gcom.xfn * uin  + gcom.yfn * vin  + gcom.zfn * win  - gcom.vfn;

    Real cref = sqrt( ABS( nscom.gama_ref * pref / rref ) );
    Real cin  = sqrt( ABS( nscom.gama    * pin  / rin  ) );

    Real gamm1 = nscom.gama - one;

    Real velin = DIST( uin, vin, win );

    //³¬ÉùËÙ
    if ( velin > cin )
    {
        if ( vnin >= 0.0 )
        {
            for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
            {
                nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
            }

            for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
            {
                nscom.primt2[ iEqu ] = 2.0 * nscom.prims1[ iEqu ] - nscom.prims2[ iEqu ];
            }

            if ( nscom.primt2[ IDX::IR ] <= 0.0 || nscom.primt2[ IDX::IP ] <= 0.0 )
            {
                for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
                {
                    nscom.primt2[ iEqu ] = nscom.prims1[ iEqu ];
                }
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
            {
                nscom.primt1[ iEqu ] = nscom.inflow[ iEqu ];
                nscom.primt2[ iEqu ] = nscom.inflow[ iEqu ];
            }
        }
    }
    else
    {
        //subsonic
        Real riemp = vnin  + 2.0 * cin  / gamm1;
        Real riemm = vnref - 2.0 * cref / gamm1;
        Real vnb   = half   * ( riemp + riemm );
        Real cb    = fourth * ( riemp - riemm ) * gamm1;

        Real vtx , vty , vtz , entr;
        if ( vnb >= 0.0 )
        {
            // exit
            entr = pin / pow( rin, nscom.gama );

            vtx = uin - gcom.xfn * vnin;
            vty = vin - gcom.yfn * vnin;
            vtz = win - gcom.zfn * vnin;
        }
        else
        {
            //inlet
            entr = pref / pow( rref, nscom.gama );
            vtx = uref - gcom.xfn * vnref;
            vty = vref - gcom.yfn * vnref;
            vtz = wref - gcom.zfn * vnref;
        }

        Real rb  = pow( ( cb * cb / ( entr * nscom.gama ) ), one / gamm1 );
        Real ub  = vtx + gcom.xfn * vnb;
        Real vb  = vty + gcom.yfn * vnb;
        Real wb  = vtz + gcom.zfn * vnb;
        Real pb  = cb * cb * rb / nscom.gama;

        nscom.primt1[ IDX::IR ] = rb;
        nscom.primt1[ IDX::IU ] = ub;
        nscom.primt1[ IDX::IV ] = vb;
        nscom.primt1[ IDX::IW ] = wb;
        nscom.primt1[ IDX::IP ] = pb;

        for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
        {
            nscom.primt2[ iEqu ] = nscom.primt1[ iEqu ];
        }
    }
}

void NsBcSolver::IsothermalVisWallBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    this->VelocityBc();

    nscom.twall = nscom.twall_dim / nscom.tref_dim;

    Real tlim = 0.1 * nscom.twall;

    Real rm = nscom.prims1[ IDX::IR ];
    Real pm = nscom.prims1[ IDX::IP ];
    Real temperature = nscom.ts1[ IDX::ITT ];

    Real rw_face = this->CalcDensity( nscom.prims1, pm, nscom.twall );
    nscom.prim[ IDX::IR ] = rw_face;
    nscom.prim[ IDX::IU ] = 0.0;
    nscom.prim[ IDX::IV ] = 0.0;
    nscom.prim[ IDX::IW ] = 0.0;
    nscom.prim[ IDX::IP ] = pm;

    Real pg1 = pm;
    Real tg1 = 2.0 * nscom.twall - temperature;

    if ( tg1 < tlim )
    {
        tg1 = tlim;
    }

    Real rg1;
    Real rw = this->CalcDensity( nscom.prims1, pg1, nscom.twall );

    rg1 = 2.0 * rw - rm;
    rg1 = MAX( 0.5 * rw, rg1 );

    nscom.primt1[ IDX::IR ] = rg1;
    nscom.primt1[ IDX::IP ] = pg1;

    nscom.primt2[ IDX::IP ] = nscom.prims2[ IDX::IP ];

    rm = nscom.prims2[ IDX::IR ];
    pm = nscom.prims2[ IDX::IP ];
    temperature = nscom.ts2[ IDX::ITT ];

    Real pg2 = pm;
    Real tg2 = 2.0 * nscom.twall - temperature;

    Real rg2;

    if ( tg2 < tlim ) tg2 = tlim;

    rg2 = this->CalcDensity( nscom.prims2, pg2, tg2 );

    nscom.primt2[ IDX::IR ] = rg2;
}

void NsBcSolver::VelocityBc()
{
    if ( nscom.bcdtkey == 0 )
    {
        nscom.primt1[ IDX::IU ] = - nscom.primt1[ IDX::IU ] + two * gcom.vfx;
        nscom.primt1[ IDX::IV ] = - nscom.primt1[ IDX::IV ] + two * gcom.vfy;
        nscom.primt1[ IDX::IW ] = - nscom.primt1[ IDX::IW ] + two * gcom.vfz;

        nscom.primt2[ IDX::IU ] = - nscom.primt2[ IDX::IU ] + two * gcom.vfx;
        nscom.primt2[ IDX::IV ] = - nscom.primt2[ IDX::IV ] + two * gcom.vfy;
        nscom.primt2[ IDX::IW ] = - nscom.primt2[ IDX::IW ] + two * gcom.vfz;
    }
    else
    {
        nscom.primt1[ IDX::IU ] = - nscom.primt1[ IDX::IU ] + two * ( * nscom.bcflow )[ IDX::IU ];
        nscom.primt1[ IDX::IV ] = - nscom.primt1[ IDX::IV ] + two * ( * nscom.bcflow )[ IDX::IV ];
        nscom.primt1[ IDX::IW ] = - nscom.primt1[ IDX::IW ] + two * ( * nscom.bcflow )[ IDX::IW ];

        nscom.primt2[ IDX::IU ] = - nscom.primt2[ IDX::IU ] + two * ( * nscom.bcflow )[ IDX::IU ];
        nscom.primt2[ IDX::IV ] = - nscom.primt2[ IDX::IV ] + two * ( * nscom.bcflow )[ IDX::IV ];
        nscom.primt2[ IDX::IW ] = - nscom.primt2[ IDX::IW ] + two * ( * nscom.bcflow )[ IDX::IW ];
    }
}

void NsBcSolver::AdiabaticVisWallBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
    this->VelocityBc();

    nscom.prim[ IDX::IR ] = nscom.prims1[ IDX::IR ];
    nscom.prim[ IDX::IU ] = 0.0;
    nscom.prim[ IDX::IV ] = 0.0;
    nscom.prim[ IDX::IW ] = 0.0;
    nscom.prim[ IDX::IP ] = nscom.prims1[ IDX::IP ];
}

void NsBcSolver::SymmetryBc()
{
    Real vx1 = nscom.prims1[ IDX::IU ];
    Real vy1 = nscom.prims1[ IDX::IV ];
    Real vz1 = nscom.prims1[ IDX::IW ];

    Real vx2 = nscom.prims2[ IDX::IU ];
    Real vy2 = nscom.prims2[ IDX::IV ];
    Real vz2 = nscom.prims2[ IDX::IW ];

    Real vnRelative1 = gcom.xfn * vx1 + gcom.yfn * vy1 + gcom.zfn * vz1 - gcom.vfn;
    Real vnRelative2 = gcom.xfn * vx2 + gcom.yfn * vy2 + gcom.zfn * vz2 - gcom.vfn;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    nscom.primt1[ IDX::IU ] = nscom.prims1[ IDX::IU ] - two * gcom.xfn * vnRelative1;
    nscom.primt1[ IDX::IV ] = nscom.prims1[ IDX::IV ] - two * gcom.yfn * vnRelative1;
    nscom.primt1[ IDX::IW ] = nscom.prims1[ IDX::IW ] - two * gcom.zfn * vnRelative1;

    nscom.primt2[ IDX::IU ] = nscom.prims2[ IDX::IU ] - two * gcom.xfn * vnRelative2;
    nscom.primt2[ IDX::IV ] = nscom.prims2[ IDX::IV ] - two * gcom.yfn * vnRelative2;
    nscom.primt2[ IDX::IW ] = nscom.prims2[ IDX::IW ] - two * gcom.zfn * vnRelative2;
}

void NsBcSolver::OversetBc()
{
}

void NsBcSolver::InterfaceBc()
{
}

void NsBcSolver::NoBc()
{
}

void NsBcSolver::UserDefinedBc()
{
}

void NsBcSolver::PeriodicBc()
{
}


Real NsBcSolver::CalcReciMolecularWeight( RealField & prim )
{
    Real reciprocalAverageMolecularWeight = one;
    return reciprocalAverageMolecularWeight;
}

Real NsBcSolver::CalcDensity( RealField & prim, Real pres, Real temperature )
{
    Real rmw = this->CalcReciMolecularWeight( prim );
    Real density = pres / ( nscom.statecoef * temperature * rmw );
    return density;
}

EndNameSpace