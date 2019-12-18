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

#include "INsUpdate.h"
#include "NsUpdate.h"
#include "UINsUpdate.h"
#include "NsInvFlux.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

INsUpdate::INsUpdate()
{
}

INsUpdate::~INsUpdate()
{
}

void INsUpdate::CmpFlowField()
{
    PrimToQ( nscom.prim0, nscom.gama, nscom.q0 );

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q[ iEqu ] = nscom.q0[ iEqu ] + nscom.dq[ iEqu ];
    }

    QToPrim( nscom.q, nscom.gama, nscom.prim, nscom.t );

    Real density  = nscom.prim[ IDX::IR ];
    Real pressure = nscom.prim[ IDX::IP ];

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        nscom.nProbe += 1;
        if ( nscom.nProbe < 2 )
        {
            this->DumpProbeInfo();
        }

        if ( ! this->WeekSolutionFix() )
        {
            this->SolutionFix();
        }
    }
}

void INsUpdate::CmpFlowFieldHyperSonic()
{
    PrimToQ( nscom.prim0, nscom.gama, nscom.q0 );

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q[ iEqu ] = nscom.q0[ iEqu ] + nscom.dq[ iEqu ];
    }

    Real density = nscom.q[ IDX::IR ];
    Real density0 = nscom.q[ IDX::IR ];
    Real dr = density - density0;
    if ( density < 0.0 )
    {
        density  = PositiveUpdate( density0, dr );
    }

    Real rd = 1.0 / density;
    Real um  = nscom.q[ IDX::IU ] * rd;
    Real vm  = nscom.q[ IDX::IV ] * rd;
    Real wm  = nscom.q[ IDX::IW ] * rd;
    Real rem = nscom.q[ IDX::IP ];
    Real v2  = ONEFLOW::SQR( um, vm, wm );

    Real reint = rem - half * density * v2;
    Real pressure = ( nscom.gama - one ) * reint;
    Real pressure0 = nscom.prim0[ IDX::IP ];

    if ( pressure < 0.0 )
    {
        Real dp = pressure - pressure0;
        pressure  = PositiveUpdate( pressure0, dp );
    }

    Real reint_new = pressure / ( nscom.gama - one );
    if ( v2 > 0.0 )
    {
        Real coef = sqrt( ABS( rem - reint_new ) * 2 / ( density * v2 ) );
        um *= coef;
        vm *= coef;
        wm *= coef;
    }

    nscom.prim[ IDX::IR ] = density;
    nscom.prim[ IDX::IU ] = um;
    nscom.prim[ IDX::IV ] = vm;
    nscom.prim[ IDX::IW ] = wm;
    nscom.prim[ IDX::IP ] = pressure;

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        nscom.nProbe += 1;
        if ( nscom.nProbe < 2 )
        {
            this->DumpProbeInfo();
        }

        if ( ! this->WeekSolutionFix() )
        {
            this->SolutionFix();
        }
    }
}

void INsUpdate::CmpFlowFieldHyperSonic_Temperature()
{
    PrimToQ( nscom.prim0, nscom.gama, nscom.q0 );

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q[ iEqu ] = nscom.q0[ iEqu ] + nscom.dq[ iEqu ];
    }

    QToPrim( nscom.q, nscom.gama, nscom.prim, nscom.t );

    Real density  = nscom.prim[ IDX::IR ];
    Real pressure = nscom.prim[ IDX::IP ];

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        nscom.nProbe += 1;
        if ( nscom.nProbe < 2 )
        {
            this->DumpProbeInfo();
        }

        if ( ! this->WeekSolutionFix() )
        {
            this->SolutionFix();
        }
    }

    if ( density > 0.0 && pressure > 0.0 )
    {
        Real rm2 = nscom.gama * SQR( nscom.mach_ref );
        Real temp_now = rm2 * pressure / density;

        Real temp_limit = 50.0;
        if ( temp_now >= temp_limit )
        {
            density = rm2 * pressure / temp_limit;
            nscom.prim[ IDX::IR ] = density;
        }

        //if ( temp_now >= temp_limit )
        //{
        //    pressure = ( density * temp_limit ) / rm2;
        //    nscom.prim[ IDX::IP ] = pressure;
        //}

        //if ( temp_now >= temp_limit )
        //{
        //    Real coef = 1.2;
        //    pressure *= coef;
        //    density = rm2 * pressure / temp_limit;
        //    nscom.prim[ IDX::IR ] = density;
        //    nscom.prim[ IDX::IP ] = pressure;
        //}

    }

}

bool INsUpdate::WeekSolutionFix()
{
    RealField fixCoefList;
    fixCoefList.push_back( 0.5 );
    fixCoefList.push_back( 0.1 );
    fixCoefList.push_back( 0.01 );
    bool flag = false;
    for ( int iFix = 0; iFix < fixCoefList.size(); ++ iFix )
    {
        Real fixCoef = fixCoefList[ iFix ];
        for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
        {
            nscom.prim[ iEqu ] += fixCoef * nscom.dq[ iEqu ];
        }

        Real density  = nscom.prim[ IDX::IR ];
        Real pressure = nscom.prim[ IDX::IP ];

        if ( ! ( density < 0.0 || pressure < 0.0 ) )
        {
            flag = true;
            break;
        }
    }
    return flag;
}

Update * CreateINsUpdate()
{
    Update * update = new UINsUpdate();
    return update;
}


EndNameSpace