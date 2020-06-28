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

#include "NsUpdate.h"
#include "UNsUpdate.h"
#include "NsInvFlux.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

NsUpdate::NsUpdate()
{
}

NsUpdate::~NsUpdate()
{
}

void NsUpdate::CalcFlowField()
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

void NsUpdate::CalcFlowFieldHyperSonic()
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

void NsUpdate::CalcFlowFieldHyperSonic_Temperature()
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

bool NsUpdate::WeekSolutionFix()
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

Update * CreateNsUpdate()
{
    Update * update = new UNsUpdate();
    return update;
}

Real PositiveUpdate( Real pn, Real dp )
{
    Real fc = 2.0;
    Real ac = 0.2;
    Real dpn = ABS( dp / pn );

    Real pn1 = pn + dp;

    if ( dpn > ac )
    {
        pn1 = pn + dp / ( 1 + fc * ( - ac + dpn ) );
    }

    return pn1;
}

void QToPrimHypersonic( RealField & q, Real gama, RealField & prim, RealField & temp )
{
    Real density = ABS( q[ IDX::IR ] );
    prim[ IDX::IR ] = density;

    Real rd = 1.0 / density;
    Real um  = q[ IDX::IU ] * rd;
    Real vm  = q[ IDX::IV ] * rd;
    Real wm  = q[ IDX::IW ] * rd;
    Real rem = q[ IDX::IP ];
    Real v2  = ONEFLOW::SQR( um, vm, wm );

    Real reint = rem - half * density * v2;

    prim[ IDX::IR ] = density;
    prim[ IDX::IU ] = um;
    prim[ IDX::IV ] = vm;
    prim[ IDX::IW ] = wm;

    Real pressure = ( gama - one ) * reint;
    prim[ IDX::IP ] = pressure;

    Real r0 = nscom.prim0[ IDX::IR ];
    Real p0 = nscom.prim0[ IDX::IP ];

    Real dr = density - r0;
    Real dp = pressure - p0;

    pressure = PositiveUpdate( p0, dp );
    density  = PositiveUpdate( r0, dr );

    //nscom.prim[ IDX::IR ] = density;
    nscom.prim[ IDX::IP ] = pressure;

}

EndNameSpace