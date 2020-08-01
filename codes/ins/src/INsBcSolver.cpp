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

#include "INsBcSolver.h"
#include "INsInvterm.h"
#include "BcData.h"
#include "INsCom.h"
#include "UCom.h"
#include "INsCtrl.h"
#include "INsIDX.h"
#include "HXMath.h"
#include "Stop.h"
#include "Boundary.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

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
    if ( vis_model.vismodel == _INVISCID )
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

void INsBcSolver::CmpFaceBc()
{
    ( this->* bcPointer )();
}

void INsBcSolver::InFlowBc()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.inflow[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.inflow[ iEqu ];
    }
}

void INsBcSolver::OutFlowBc()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }
}

void INsBcSolver::PoleBc()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }
}

void INsBcSolver::FarFieldBc()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }

    //inner point
    Real rin, uin, vin, win, pin;
    ONEFLOW::INsExtract( inscom.prims1, rin, uin, vin, win, pin );

    gcom.xfn *= inscom.faceOuterNormal;
    gcom.yfn *= inscom.faceOuterNormal;
    gcom.zfn *= inscom.faceOuterNormal;

    Real rref = inscom.inflow[ IIDX::IIR ];
    Real uref = inscom.inflow[ IIDX::IIU ];
    Real vref = inscom.inflow[ IIDX::IIV ];
    Real wref = inscom.inflow[ IIDX::IIW ];
    Real pref = inscom.inflow[ IIDX::IIP ];

    Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
    Real vnin  = gcom.xfn * uin  + gcom.yfn * vin  + gcom.zfn * win  - gcom.vfn;

    Real cref = sqrt( ABS( inscom.gama_ref * pref / rref ) );
    Real cin  = sqrt( ABS( inscom.gama    * pin  / rin  ) );

    Real gamm1 = inscom.gama - one;

    Real velin = DIST( uin, vin, win );

    //³¬ÉùËÙ
    if ( velin > cin )
    {
        if ( vnin >= 0.0 )
        {
            for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
            {
                inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
            }

            for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
            {
                inscom.primt2[ iEqu ] = 2.0 * inscom.prims1[ iEqu ] - inscom.prims2[ iEqu ];
            }

            if ( inscom.primt2[ IIDX::IIR ] <= 0.0 || inscom.primt2[ IIDX::IIP ] <= 0.0 )
            {
                for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
                {
                    inscom.primt2[ iEqu ] = inscom.prims1[ iEqu ];
                }
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
            {
                inscom.primt1[ iEqu ] = inscom.inflow[ iEqu ];
                inscom.primt2[ iEqu ] = inscom.inflow[ iEqu ];
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
            entr = pin / pow( rin, inscom.gama );

            vtx = uin - gcom.xfn * vnin;
            vty = vin - gcom.yfn * vnin;
            vtz = win - gcom.zfn * vnin;
        }
        else
        {
            //inlet
            entr = pref / pow( rref, inscom.gama );
            vtx = uref - gcom.xfn * vnref;
            vty = vref - gcom.yfn * vnref;
            vtz = wref - gcom.zfn * vnref;
        }

        Real rb  = pow( ( cb * cb / ( entr * inscom.gama ) ), one / gamm1 );
        Real ub  = vtx + gcom.xfn * vnb;
        Real vb  = vty + gcom.yfn * vnb;
        Real wb  = vtz + gcom.zfn * vnb;
        Real pb  = cb * cb * rb / inscom.gama;

        inscom.primt1[ IIDX::IIR ] = rb;
        inscom.primt1[ IIDX::IIU ] = ub;
        inscom.primt1[ IIDX::IIV ] = vb;
        inscom.primt1[ IIDX::IIW ] = wb;
        inscom.primt1[ IIDX::IIP ] = pb;

        for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
        {
            inscom.primt2[ iEqu ] = inscom.primt1[ iEqu ];
        }
    }
}

void INsBcSolver::IsothermalVisWallBc()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }

    this->VelocityBc();

    inscom.twall = inscom.twall_dim / inscom.tref_dim;

    Real tlim = 0.1 * inscom.twall;

    Real rm = inscom.prims1[ IIDX::IIR ];
    Real pm = inscom.prims1[ IIDX::IIP ];
    Real temperature = inscom.ts1[ IIDX::IITT ];

    Real rw_face = this->CmpDensity( inscom.prims1, pm, inscom.twall );
    inscom.prim[ IIDX::IIR ] = rw_face;
    inscom.prim[ IIDX::IIU ] = 0.0;
    inscom.prim[ IIDX::IIV ] = 0.0;
    inscom.prim[ IIDX::IIW ] = 0.0;
    inscom.prim[ IIDX::IIP ] = pm;

    Real pg1 = pm;
    Real tg1 = 2.0 * inscom.twall - temperature;

    if ( tg1 < tlim )
    {
        tg1 = tlim;
    }

    Real rg1;
    Real rw = this->CmpDensity( inscom.prims1, pg1, inscom.twall );

    rg1 = 2.0 * rw - rm;
    rg1 = MAX( 0.5 * rw, rg1 );

    inscom.primt1[ IIDX::IIR ] = rg1;
    inscom.primt1[ IIDX::IIP ] = pg1;

    inscom.primt2[ IIDX::IIP ] = inscom.prims2[ IIDX::IIP ];

    rm = inscom.prims2[ IIDX::IIR ];
    pm = inscom.prims2[ IIDX::IIP ];
    temperature = inscom.ts2[ IIDX::IITT ];

    Real pg2 = pm;
    Real tg2 = 2.0 * inscom.twall - temperature;

    Real rg2;

    if ( tg2 < tlim ) tg2 = tlim;

    rg2 = this->CmpDensity( inscom.prims2, pg2, tg2 );

    inscom.primt2[ IIDX::IIR ] = rg2;
}

void INsBcSolver::VelocityBc()
{
	
    if ( inscom.bcdtkey == 0 )
    {
        inscom.primt1[ IIDX::IIU ] = - inscom.primt1[ IIDX::IIU ] + two * gcom.vfx;
        inscom.primt1[ IIDX::IIV ] = - inscom.primt1[ IIDX::IIV ] + two * gcom.vfy;
        inscom.primt1[ IIDX::IIW ] = - inscom.primt1[ IIDX::IIW ] + two * gcom.vfz;

        inscom.primt2[ IIDX::IIU ] = - inscom.primt2[ IIDX::IIU ] + two * gcom.vfx;
        inscom.primt2[ IIDX::IIV ] = - inscom.primt2[ IIDX::IIV ] + two * gcom.vfy;
        inscom.primt2[ IIDX::IIW ] = - inscom.primt2[ IIDX::IIW ] + two * gcom.vfz;
    }
    else
    {
        inscom.primt1[ IIDX::IIU ] = - inscom.primt1[ IIDX::IIU ] + two * ( * inscom.bcflow )[ IIDX::IIU ];
        inscom.primt1[ IIDX::IIV ] = - inscom.primt1[ IIDX::IIV ] + two * ( * inscom.bcflow )[ IIDX::IIV ];
        inscom.primt1[ IIDX::IIW ] = - inscom.primt1[ IIDX::IIW ] + two * ( * inscom.bcflow )[ IIDX::IIW ];

        inscom.primt2[ IIDX::IIU ] = - inscom.primt2[ IIDX::IIU ] + two * ( * inscom.bcflow )[ IIDX::IIU ];
        inscom.primt2[ IIDX::IIV ] = - inscom.primt2[ IIDX::IIV ] + two * ( * inscom.bcflow )[ IIDX::IIV ];
        inscom.primt2[ IIDX::IIW ] = - inscom.primt2[ IIDX::IIW ] + two * ( * inscom.bcflow )[ IIDX::IIW ];
    }
}

void INsBcSolver::AdiabaticVisWallBc()
{

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }
    this->VelocityBc();

    inscom.prim[ IIDX::IIR ] = inscom.prims1[ IIDX::IIR ];
    inscom.prim[ IIDX::IIU ] = 0.0;
    inscom.prim[ IIDX::IIV ] = 0.0;
    inscom.prim[ IIDX::IIW ] = 0.0;
    inscom.prim[ IIDX::IIP ] = inscom.prims1[ IIDX::IIP ];
}

void INsBcSolver::SymmetryBc()
{
    Real vx1 = inscom.prims1[ IIDX::IIU ];
    Real vy1 = inscom.prims1[ IIDX::IIV ];
    Real vz1 = inscom.prims1[ IIDX::IIW ];

    Real vx2 = inscom.prims2[ IIDX::IIU ];
    Real vy2 = inscom.prims2[ IIDX::IIV ];
    Real vz2 = inscom.prims2[ IIDX::IIW ];

    Real vnRelative1 = gcom.xfn * vx1 + gcom.yfn * vy1 + gcom.zfn * vz1 - gcom.vfn;
    Real vnRelative2 = gcom.xfn * vx2 + gcom.yfn * vy2 + gcom.zfn * vz2 - gcom.vfn;

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.primt1[ iEqu ] = inscom.prims1[ iEqu ];
        inscom.primt2[ iEqu ] = inscom.prims2[ iEqu ];
    }

    inscom.primt1[ IIDX::IIU ] = inscom.prims1[ IIDX::IIU ] - two * gcom.xfn * vnRelative1;
    inscom.primt1[ IIDX::IIV ] = inscom.prims1[ IIDX::IIV ] - two * gcom.yfn * vnRelative1;
    inscom.primt1[ IIDX::IIW ] = inscom.prims1[ IIDX::IIW ] - two * gcom.zfn * vnRelative1;

    inscom.primt2[ IIDX::IIU ] = inscom.prims2[ IIDX::IIU ] - two * gcom.xfn * vnRelative2;
    inscom.primt2[ IIDX::IIV ] = inscom.prims2[ IIDX::IIV ] - two * gcom.yfn * vnRelative2;
    inscom.primt2[ IIDX::IIW ] = inscom.prims2[ IIDX::IIW ] - two * gcom.zfn * vnRelative2;
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


Real INsBcSolver::CmpReciMolecularWeight( RealField & prim )
{
    Real reciprocalAverageMolecularWeight = one;
    return reciprocalAverageMolecularWeight;
}

Real INsBcSolver::CmpDensity( RealField & prim, Real pres, Real temperature )
{
    Real rmw = this->CmpReciMolecularWeight( prim );
    Real density = pres / ( inscom.statecoef * temperature * rmw );
    return density;
}

EndNameSpace

