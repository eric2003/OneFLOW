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

#include "INsCom.h"
#include "INsIdx.h"
#include "DataBase.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Chemical.h"

BeginNameSpace( ONEFLOW )

INsCom inscom;

INsCom::INsCom()
{
    init_flag = false;
}

INsCom::~INsCom()
{
    ;
}

void INsCom::Init()
{
    if ( init_flag ) return;
    init_flag = true;
    icmpInv = 1;
    ischeme = GetDataValue< int >( "ischeme" );
    ivischeme = GetDataValue< int >( "ivischeme" );
    timestepModel = GetDataValue< int >( "timestepModel" );
    visSRModel = GetDataValue< int >( "visSRModel" );
    visTimestepModel = GetDataValue< int >( "visTimestepModel" );
    chemModel = GetDataValue< int >( "chemModel" );
    nTModel = GetDataValue< int >( "nTModel" );

    faceOuterNormal = 1;
    reynolds = GetDataValue< Real >( "reynolds" );

    prl = GetDataValue< Real >( "prl" );
    prt = GetDataValue< Real >( "prt" );

    aoa = GetDataValue< Real >( "aoa_degree" ) * PI / 180.0;
    aos = GetDataValue< Real >( "sideslip_degree" ) * PI / 180.0;

    mach_ref = GetDataValue< Real >( "mach_ref" );
    gama_ref = GetDataValue< Real >( "gama_ref" );

    twall_dim = GetDataValue< Real >( "twall_dim" );
    tref_dim = GetDataValue< Real >( "tref_dim" );
    pref_dim = GetDataValue< Real >( "pref_dim" );
    dref_dim = GetDataValue< Real >( "dref_dim" );
    vref_dim = GetDataValue< Real >( "vref_dim" );

    elevation = GetDataValue< Real >( "elevation" );
    reylref_dim = GetDataValue< Real >( "reylref_dim" );

    gasInfoStrategy = GetDataValue< int >( "gasInfoStrategy" );
    gasModelFile = GetDataValue< string >( "gasModelFile" );
    machStrategy = GetDataValue< int >( "machStrategy" );

    schmidtl = GetDataValue< Real >( "schmidtl" );
    schmidtt = GetDataValue< Real >( "schmidtt" );

    max_time_ratio = GetDataValue< Real >( "max_time_ratio" );

    nEqu = GetDataValue< int >( "nEqu" );
    nTEqu = GetDataValue< int >( "nTEqu" );
    chem.INsInit();

    oprl = 1.0 / prl;
    oprt = 1.0 / prt;

    twall = twall_dim / tref_dim;

    const_cp = 1.0 / ( ( gama_ref - 1.0 ) * SQR( mach_ref ) );

    q1.resize( nTEqu );
    q2.resize( nTEqu );

    q.resize( nTEqu );
    q0.resize( nTEqu );
    dq.resize( nTEqu );

    prim.resize( nTEqu );
    prim0.resize( nTEqu );

    prims1.resize( nTEqu );
    prims2.resize( nTEqu );

    primt1.resize( nTEqu );
    primt2.resize( nTEqu );

    t.resize( nTModel );
    t0.resize( nTModel );

    ts1.resize( nTModel );
    ts2.resize( nTModel );

    tt1.resize( nTModel );
    tt2.resize( nTModel );
}

void INsExtract(RealField & prim, Real & rm, Real & um, Real & vm, Real & wm, Real & pm)
{
	rm = prim[IIDX::IIR];
	um = prim[IIDX::IIU];
	vm = prim[IIDX::IIV];
	wm = prim[IIDX::IIW];
	pm = prim[IIDX::IIP];
}


bool INsCheckFunction( RealField & q )
{
    if ( q[ IIDX::IIR ] <= 0.0 || q[ IIDX::IIP ] <= 0.0 ) return false;
    return true;
}

EndNameSpace