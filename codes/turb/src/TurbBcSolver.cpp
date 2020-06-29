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

#include "TurbBcSolver.h"
#include "BcData.h"
#include "UCom.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "Ctrl.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Stop.h"
#include "Boundary.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

TurbBcSolver::TurbBcSolver()
{
    ;
}

TurbBcSolver::~TurbBcSolver()
{
    ;
}

void TurbBcSolver::SetBc()
{
    this->bcPointer = & TurbBcSolver::NoBc;

    this->updateFlag = true;

    if ( ug.bctype < 0 )
    {
        this->bcPointer  = & TurbBcSolver::InterfaceBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::OVERSET )
    {
        this->bcPointer  = & TurbBcSolver::OversetBc;
        this->updateFlag = false;
    }
    else if ( ug.bctype == BC::EXTRAPOLATION )
    {
        this->bcPointer = & TurbBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::SOLID_SURFACE )
    {
        this->bcPointer = & TurbBcSolver::VisWallBc;
    }
    else if ( ug.bctype == BC::SYMMETRY )
    {
        this->bcPointer = & TurbBcSolver::SymmetryBc;
    }
    else if ( ug.bctype == BC::FARFIELD )
    {
        this->bcPointer = & TurbBcSolver::FarFieldBc;
    }
    else if ( ug.bctype == BC::INFLOW )
    {
        this->bcPointer = & TurbBcSolver::InFlowBc;
    }
    else if ( ug.bctype == BC::OUTFLOW )
    {
        this->bcPointer = & TurbBcSolver::OutFlowBc;
    }
    else if ( ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE )
    {
        this->bcPointer = & TurbBcSolver::PoleBc;
    }
    else if ( ug.bctype == BC::GENERIC_2 )
    {
        //userDefined
        this->bcPointer = & TurbBcSolver::UserDefinedBc;
    }
    else if ( ug.bctype == BC::PERIODIC ) 
    {
        this->bcPointer = & TurbBcSolver::PeriodicBc;
    }
    else
    {
        cout << "Error : Illegal BCtype ID " << ug.bctype << endl;
        Stop("");
    }
}

void TurbBcSolver::CalcFaceBc()
{
    ( this->* bcPointer )();
}

void TurbBcSolver::InFlowBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.inflow[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.inflow[ iEqu ];
    }
}

void TurbBcSolver::OutFlowBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.prims2[ iEqu ];
    }
}

void TurbBcSolver::PoleBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.prims2[ iEqu ];
    }
}

void TurbBcSolver::FarFieldBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.prims2[ iEqu ];
    }

    //inner point
    Real rin, uin, vin, win, pin;
    ONEFLOW::Extract( turbcom.ns_prims1, rin, uin, vin, win, pin );

    gcom.xfn *= gcom.faceOuterNormal;
    gcom.yfn *= gcom.faceOuterNormal;
    gcom.zfn *= gcom.faceOuterNormal;

    Real rref = nscom.inflow[ IDX::IR ];
    Real uref = nscom.inflow[ IDX::IU ];
    Real vref = nscom.inflow[ IDX::IV ];
    Real wref = nscom.inflow[ IDX::IW ];
    Real pref = nscom.inflow[ IDX::IP ];

    Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
    Real vnin  = gcom.xfn * uin  + gcom.yfn * vin  + gcom.zfn * win  - gcom.vfn;

    Real cref = sqrt( ABS( nscom.gama_ref * pref / rref ) );
    Real cin  = sqrt( ABS( nscom.gama    * pin  / rin  ) );

    Real velin = DIST( uin, vin, win );

    //超声速
    if ( velin > cin )
    {
        if ( vnin >= 0.0 )
        {
            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
                turbcom.primt2[ iEqu ] = turbcom.prims1[ iEqu ];
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                turbcom.primt1[ iEqu ] = turbcom.inflow[ iEqu ];
                turbcom.primt2[ iEqu ] = turbcom.inflow[ iEqu ];
            }
        }
    }
    else
    {
        Real gamm1 = nscom.gama - one;
        //subsonic
        Real riemp = vnin  + 2.0 * cin  / gamm1;
        Real riemm = vnref - 2.0 * cref / gamm1;
        Real vnb   = half   * ( riemp + riemm );
        Real cb    = fourth * ( riemp - riemm ) * gamm1;

        if ( vnb >= 0.0 )
        {
            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
                turbcom.primt2[ iEqu ] = turbcom.prims1[ iEqu ];
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                turbcom.primt1[ iEqu ] = turbcom.inflow[ iEqu ];
                turbcom.primt2[ iEqu ] = turbcom.inflow[ iEqu ];
            }
        }
    }
}

void TurbBcSolver::VisWallBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.prims2[ iEqu ];
    }

    if ( turbcom.nEqu == 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turbcom.primt1[ iEqu ] = - turbcom.prims1[ iEqu ];
            turbcom.primt2[ iEqu ] = - turbcom.prims2[ iEqu ];
        }
    }
    else if ( turbcom.nEqu >= 2 )  // 考虑转捩方程后，多了两方程
    {
        Real dist2 = SQR( turbcom.dist );

        Real kwWall = 60.0 * turbcom.visl / ( turbcom.rho * turbcom.beta1 * dist2 * turbcom.reynolds );

        turbcom.primt1[ IKE ] = - turbcom.prims1[ IKE ];
        turbcom.primt2[ IKE ] = - turbcom.prims2[ IKE ];

        //turbcom.primt1[ IKW ] = MAX( two * kwWall - turbcom.prims1[ IKW ], SMALL );
        //turbcom.primt2[ IKW ] = MAX( two * kwWall - turbcom.prims2[ IKW ], SMALL );

        turbcom.primt1[ IKW ] = kwWall;
        turbcom.primt2[ IKW ] = kwWall;
    }

}

void TurbBcSolver::SymmetryBc()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.primt1[ iEqu ] = turbcom.prims1[ iEqu ];
        turbcom.primt2[ iEqu ] = turbcom.prims2[ iEqu ];
    }
}

void TurbBcSolver::OversetBc()
{
}

void TurbBcSolver::InterfaceBc()
{
}

void TurbBcSolver::NoBc()
{
}

void TurbBcSolver::UserDefinedBc()
{
}

void TurbBcSolver::PeriodicBc()
{
}

EndNameSpace