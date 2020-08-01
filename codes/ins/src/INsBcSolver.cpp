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

#include "INsBcSolver.h"
#include "BcData.h"
#include "FlowModel.h"
#include "INsCom.h"
#include "UCom.h"
#include "INsCtrl.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Stop.h"
#include "Boundary.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

BcData ins_bc_data;

INsBcSolver::INsBcSolver()
{
    ;
}

INsBcSolver::~INsBcSolver()
{
    ;
}

void INsBcSolver::SetBc()
{
    this->bcPointer = & INsBcSolver::NoBc;

    this->updateFlag = true;

    if ( ug.bctype < 0 )
    {
        this->bcPointer  = & INsBcSolver::InterfaceBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::OVERSET )
    {
        this->bcPointer  = & INsBcSolver::OversetBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::EXTRAPOLATION )
    {
        this->bcPointer = & INsBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::SOLID_SURFACE )
    {
        this->SetSolidSurfaceBc();
    }
    else if ( ug.bctype == BC::SYMMETRY )
    {
        this->bcPointer = & INsBcSolver::SymmetryBc;
    }
    else if ( ug.bctype == BC::FARFIELD )
    {
        this->bcPointer = & INsBcSolver::FarFieldBc;
    }
    else if ( ug.bctype == BC::INFLOW )
    {
        this->bcPointer = & INsBcSolver::InFlowBc;
    }
    else if ( ug.bctype == BC::OUTFLOW )
    {
        this->bcPointer = & INsBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE )
    {
        this->bcPointer = & INsBcSolver::PoleBc;
    }
    else if ( ug.bctype == BC::GENERIC_2 )
    {
        //userDefined
        this->bcPointer = & INsBcSolver::UserDefinedBc;
    }
    else if ( ug.bctype == BC::PERIODIC ) 
    {
        this->bcPointer = & INsBcSolver::PeriodicBc;
    }
    else
    {
        cout << "Error : Illegal BCtype ID " << ug.bctype << endl;
        Stop("");
    }
}

void INsBcSolver::SetSolidSurfaceBc()
{
    if ( vis_model.vismodel == INVISCID )
    {
        this->bcPointer = & INsBcSolver::SymmetryBc;
    }
    else if ( ins_ctrl.isowallbc == 0 )
    {
        //viscous adiabatic wall
        this->bcPointer = & INsBcSolver::AdiabaticVisWallBc;
    }
    else
    {
        //viscous iso-thermal wall
        this->bcPointer = & INsBcSolver::IsothermalVisWallBc;
    }
}

void INsBcSolver::CalcFaceBc()
{
    ( this->* bcPointer )();
}

void INsBcSolver::InFlowBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.inflow[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.inflow[ iEqu ];
    }
}

void INsBcSolver::OutFlowBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
}

void INsBcSolver::PoleBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
}

void INsBcSolver::FarFieldBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    //inner point
    Real rin, uin, vin, win, pin;
    ONEFLOW::INsExtract( nscom.prims1, rin, uin, vin, win, pin );

    gcom.xfn *= nscom.faceOuterNormal;
    gcom.yfn *= nscom.faceOuterNormal;
    gcom.zfn *= nscom.faceOuterNormal;

    Real rref = nscom.inflow[ IIDX::IIR ];
    Real uref = nscom.inflow[ IIDX::IIU ];
    Real vref = nscom.inflow[ IIDX::IIV ];
    Real wref = nscom.inflow[ IIDX::IIW ];
    Real pref = nscom.inflow[ IIDX::IIP ];

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

            if ( nscom.primt2[ IIDX::IIR ] <= 0.0 || nscom.primt2[ IIDX::IIP ] <= 0.0 )
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

        nscom.primt1[ IIDX::IIR ] = rb;
        nscom.primt1[ IIDX::IIU ] = ub;
        nscom.primt1[ IIDX::IIV ] = vb;
        nscom.primt1[ IIDX::IIW ] = wb;
        nscom.primt1[ IIDX::IIP ] = pb;

        for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
        {
            nscom.primt2[ iEqu ] = nscom.primt1[ iEqu ];
        }
    }
}

void INsBcSolver::IsothermalVisWallBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    this->VelocityBc();

    nscom.twall = nscom.twall_dim / nscom.tref_dim;

    Real tlim = 0.1 * nscom.twall;

    Real rm = nscom.prims1[ IIDX::IIR ];
    Real pm = nscom.prims1[ IIDX::IIP ];
    Real temperature = nscom.ts1[ IIDX::IITT ];

    Real rw_face = this->CalcDensity( nscom.prims1, pm, nscom.twall );
    nscom.prim[ IIDX::IIR ] = rw_face;
    nscom.prim[ IIDX::IIU ] = 0.0;
    nscom.prim[ IIDX::IIV ] = 0.0;
    nscom.prim[ IIDX::IIW ] = 0.0;
    nscom.prim[ IIDX::IIP ] = pm;

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

    nscom.primt1[ IIDX::IIR ] = rg1;
    nscom.primt1[ IIDX::IIP ] = pg1;

    nscom.primt2[ IIDX::IIP ] = nscom.prims2[ IIDX::IIP ];

    rm = nscom.prims2[ IIDX::IIR ];
    pm = nscom.prims2[ IIDX::IIP ];
    temperature = nscom.ts2[ IIDX::IITT ];

    Real pg2 = pm;
    Real tg2 = 2.0 * nscom.twall - temperature;

    Real rg2;

    if ( tg2 < tlim ) tg2 = tlim;

    rg2 = this->CalcDensity( nscom.prims2, pg2, tg2 );

    nscom.primt2[ IIDX::IIR ] = rg2;
}

void INsBcSolver::VelocityBc()
{
    if ( nscom.bcdtkey == 0 )
    {
        nscom.primt1[ IIDX::IIU ] = - nscom.primt1[ IIDX::IIU ] + two * gcom.vfx;
        nscom.primt1[ IIDX::IIV ] = - nscom.primt1[ IIDX::IIV ] + two * gcom.vfy;
        nscom.primt1[ IIDX::IIW ] = - nscom.primt1[ IIDX::IIW ] + two * gcom.vfz;

        nscom.primt2[ IIDX::IIU ] = - nscom.primt2[ IIDX::IIU ] + two * gcom.vfx;
        nscom.primt2[ IIDX::IIV ] = - nscom.primt2[ IIDX::IIV ] + two * gcom.vfy;
        nscom.primt2[ IIDX::IIW ] = - nscom.primt2[ IIDX::IIW ] + two * gcom.vfz;
    }
    else
    {
        nscom.primt1[ IIDX::IIU ] = - nscom.primt1[ IIDX::IIU ] + two * ( * nscom.bcflow )[ IIDX::IIU ];
        nscom.primt1[ IIDX::IIV ] = - nscom.primt1[ IIDX::IIV ] + two * ( * nscom.bcflow )[ IIDX::IIV ];
        nscom.primt1[ IIDX::IIW ] = - nscom.primt1[ IIDX::IIW ] + two * ( * nscom.bcflow )[ IIDX::IIW ];

        nscom.primt2[ IIDX::IIU ] = - nscom.primt2[ IIDX::IIU ] + two * ( * nscom.bcflow )[ IIDX::IIU ];
        nscom.primt2[ IIDX::IIV ] = - nscom.primt2[ IIDX::IIV ] + two * ( * nscom.bcflow )[ IIDX::IIV ];
        nscom.primt2[ IIDX::IIW ] = - nscom.primt2[ IIDX::IIW ] + two * ( * nscom.bcflow )[ IIDX::IIW ];
    }
}

void INsBcSolver::AdiabaticVisWallBc()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }
    this->VelocityBc();

    nscom.prim[ IIDX::IIR ] = nscom.prims1[ IIDX::IIR ];
    nscom.prim[ IIDX::IIU ] = 0.0;
    nscom.prim[ IIDX::IIV ] = 0.0;
    nscom.prim[ IIDX::IIW ] = 0.0;
    nscom.prim[ IIDX::IIP ] = nscom.prims1[ IIDX::IIP ];
}

void INsBcSolver::SymmetryBc()
{
    Real vx1 = nscom.prims1[ IIDX::IIU ];
    Real vy1 = nscom.prims1[ IIDX::IIV ];
    Real vz1 = nscom.prims1[ IIDX::IIW ];

    Real vx2 = nscom.prims2[ IIDX::IIU ];
    Real vy2 = nscom.prims2[ IIDX::IIV ];
    Real vz2 = nscom.prims2[ IIDX::IIW ];

    Real vnRelative1 = gcom.xfn * vx1 + gcom.yfn * vy1 + gcom.zfn * vz1 - gcom.vfn;
    Real vnRelative2 = gcom.xfn * vx2 + gcom.yfn * vy2 + gcom.zfn * vz2 - gcom.vfn;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.primt1[ iEqu ] = nscom.prims1[ iEqu ];
        nscom.primt2[ iEqu ] = nscom.prims2[ iEqu ];
    }

    nscom.primt1[ IIDX::IIU ] = nscom.prims1[ IIDX::IIU ] - two * gcom.xfn * vnRelative1;
    nscom.primt1[ IIDX::IIV ] = nscom.prims1[ IIDX::IIV ] - two * gcom.yfn * vnRelative1;
    nscom.primt1[ IIDX::IIW ] = nscom.prims1[ IIDX::IIW ] - two * gcom.zfn * vnRelative1;

    nscom.primt2[ IIDX::IIU ] = nscom.prims2[ IIDX::IIU ] - two * gcom.xfn * vnRelative2;
    nscom.primt2[ IIDX::IIV ] = nscom.prims2[ IIDX::IIV ] - two * gcom.yfn * vnRelative2;
    nscom.primt2[ IIDX::IIW ] = nscom.prims2[ IIDX::IIW ] - two * gcom.zfn * vnRelative2;
}

void INsBcSolver::OversetBc()
{
}

void INsBcSolver::InterfaceBc()
{
}

void INsBcSolver::NoBc()
{
}

void INsBcSolver::UserDefinedBc()
{
}

void INsBcSolver::PeriodicBc()
{
}


Real INsBcSolver::CalcReciMolecularWeight( RealField & prim )
{
    Real reciprocalAverageMolecularWeight = one;
    return reciprocalAverageMolecularWeight;
}

Real INsBcSolver::CalcDensity( RealField & prim, Real pres, Real temperature )
{
    Real rmw = this->CalcReciMolecularWeight( prim );
    Real density = pres / ( nscom.statecoef * temperature * rmw );
    return density;
}

EndNameSpace

