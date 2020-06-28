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

#include "TurbLusgs.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include "UCom.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

TurbLusgsData turblu;

TurbLusgsData::TurbLusgsData()
{
    ;
}

TurbLusgsData::~TurbLusgsData()
{
    ;
}

void TurbLusgsData::Init()
{
    nEqu = turbcom.nEqu;
    radius.resize( nEqu );
    matrix.resize( nEqu );
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

TurbLusgs::TurbLusgs()
{
}

TurbLusgs::~TurbLusgs()
{
}

void TurbLusgs::InitializeSub()
{
}

void TurbLusgs::DumpSweepInformation()
{
    int pid = ONEFLOW::Parallel::GetPid();
    if ( pid == ONEFLOW::Parallel::GetServerid() )
    {
        cout << "Turb residual reduced by " << turblu.dmax << " with " << turblu.numberOfRealSweeps << " Sweeps\n";
    }
}

void TurbLusgs::ZeroFluxIncrement()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.rhs0[ iEqu ] = 0.0;
    }
}

void TurbLusgs::AddFluxIncrement()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.rhs0[ iEqu ] += turblu.dfj[ iEqu ];
    }

    //if ( ug.cId == 22 )
    //{
    //    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    //    {
    //        cout << " add " << turblu.rhs0[ iEqu ] << " " << turblu.dfj[ iEqu ] << "\n";
    //    }
    //}
    
}

void TurbLusgs::AddFluxIncrement( const Real & coef )
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.rhs0[ iEqu ] += coef * turblu.dfj[ iEqu ];
    }
}

void TurbLusgs::GetFluxIncrement( int signOfMatrix )
{
    this->GetStandardFluxIncrement( signOfMatrix );
}

void TurbLusgs::GetStandardFluxIncrement( int signOfMatrix )
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.dfj[ iEqu ] = turblu.matrix[ iEqu ] * turblu.dqj[ iEqu ];
    }
}

void TurbLusgs::InitializeSweep( int iSweep )
{
    turblu.norm = 0.0;
}

bool TurbLusgs::UpdateSweep( int iSweep )
{
    turblu.numberOfRealSweeps = iSweep + 1;
    if ( iSweep == 0 )
    {
        turblu.norm0 = turblu.norm;
        turblu.dmax  = 1.0;
    }
    else
    {
        turblu.dmax = sqrt( turblu.norm / ( turblu.norm0 + SMALL ) );
    }

    if ( turblu.dmax < turblu.tol )
    {
        return true;
    }
    return false;
}

void TurbLusgs::CalcLowerChange()
{
    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] = turblu.dqi[ iEqu ] - turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] /=  turblu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi[ iEqu ] = turblu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.drhs[ iEqu ] += turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqSweep =  turblu.dqi[ iEqu ] - turblu.dqi0[ iEqu ];
            turblu.norm   += SQR( turblu.dqSweep );
        }
    }
    else
    {

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] = turblu.rhs[ iEqu ] - turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] /=  turblu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi[ iEqu ] = turblu.tmp[ iEqu ];
        }
    }
}

void TurbLusgs::CalcUpperChange()
{
    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] = turblu.dqi[ iEqu ] - turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] /=  turblu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi[ iEqu ] = turblu.tmp[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {            
            turblu.drhs[ iEqu ] += turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqSweep     = turblu.dqi[ iEqu ] - turblu.dqi0[ iEqu ];
            turblu.norm       += SQR( turblu.dqSweep );
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] = - turblu.rhs0[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.tmp[ iEqu ] /=  turblu.radius[ iEqu ];
        }

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi[ iEqu ] += turblu.tmp[ iEqu ];
        }
    }
}

bool TurbLusgs::IsOversetCell()
{
    return ( gcom.blank <= 0 );
}

void TurbLusgs::ZeroOversetCell()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.dqi[ iEqu ] = 0.0;
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.drhs[ iEqu ] = 0.0;
    }
}


EndNameSpace