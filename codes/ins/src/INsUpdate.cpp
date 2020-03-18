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
#include "INsInvterm.h"
#include "INsCom.h"
#include "INsIDX.h"
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
    //INsPrimToQ( inscom.prim0, inscom.gama, inscom.q0 );

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.q[ iEqu ] = inscom.q0[ iEqu ] + inscom.dq[ iEqu ];
    }

    //INsQToPrim( inscom.q, inscom.gama, inscom.prim, inscom.t );

    Real density  = inscom.prim[ IIDX::IIR ];
    Real pressure = inscom.prim[ IIDX::IIP ];

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        inscom.nProbe += 1;
        if ( inscom.nProbe < 2 )
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
   // INsPrimToQ( inscom.prim0, inscom.gama, inscom.q0 );

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.q[ iEqu ] = inscom.q0[ iEqu ] + inscom.dq[ iEqu ];
    }

    Real density = inscom.q[ IIDX::IIR ];
    Real density0 = inscom.q[ IIDX::IIR ];
    Real dr = density - density0;
    if ( density < 0.0 )
    {
        density  = PositiveUpdate( density0, dr );
    }

    Real rd = 1.0 / density;
    Real um  = inscom.q[ IIDX::IIU ] * rd;
    Real vm  = inscom.q[ IIDX::IIV ] * rd;
    Real wm  = inscom.q[ IIDX::IIW ] * rd;
    Real rem = inscom.q[ IIDX::IIP ];
    Real v2  = ONEFLOW::SQR( um, vm, wm );

    Real reint = rem - half * density * v2;
    Real pressure = ( inscom.gama - one ) * reint;
    Real pressure0 = inscom.prim0[ IIDX::IIP ];

    if ( pressure < 0.0 )
    {
        Real dp = pressure - pressure0;
        pressure  = PositiveUpdate( pressure0, dp );
    }

    Real reint_new = pressure / ( inscom.gama - one );
    if ( v2 > 0.0 )
    {
        Real coef = sqrt( ABS( rem - reint_new ) * 2 / ( density * v2 ) );
        um *= coef;
        vm *= coef;
        wm *= coef;
    }

    inscom.prim[ IIDX::IIR ] = density;
    inscom.prim[ IIDX::IIU ] = um;
    inscom.prim[ IIDX::IIV ] = vm;
    inscom.prim[ IIDX::IIW ] = wm;
    inscom.prim[ IIDX::IIP ] = pressure;

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        inscom.nProbe += 1;
        if ( inscom.nProbe < 2 )
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
//	INsPrimToQ( inscom.prim0, inscom.gama, inscom.q0 );

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.q[ iEqu ] = inscom.q0[ iEqu ] + inscom.dq[ iEqu ];
    }

//	INsQToPrim( inscom.q, inscom.gama, inscom.prim, inscom.t );

    Real density  = inscom.prim[ IIDX::IIR ];
    Real pressure = inscom.prim[ IIDX::IIP ];

    if ( density <= 0.0 || pressure <= 0.0 )
    {
        inscom.nProbe += 1;
        if ( inscom.nProbe < 2 )
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
        Real rm2 = inscom.gama * SQR( inscom.mach_ref );
        Real temp_now = rm2 * pressure / density;

        Real temp_limit = 50.0;
        if ( temp_now >= temp_limit )
        {
            density = rm2 * pressure / temp_limit;
            inscom.prim[ IIDX::IIR ] = density;
        }

        //if ( temp_now >= temp_limit )
        //{
        //    pressure = ( density * temp_limit ) / rm2;
        //    nscom.prim[ IIDX::IP ] = pressure;
        //}

        //if ( temp_now >= temp_limit )
        //{
        //    Real coef = 1.2;
        //    pressure *= coef;
        //    density = rm2 * pressure / temp_limit;
        //    nscom.prim[ IIDX::IR ] = density;
        //    nscom.prim[ IIDX::IP ] = pressure;
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
        for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
        {
            inscom.prim[ iEqu ] += fixCoef * inscom.dq[ iEqu ];
        }

        Real density  = inscom.prim[ IIDX::IIR ];
        Real pressure = inscom.prim[ IIDX::IIP ];

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