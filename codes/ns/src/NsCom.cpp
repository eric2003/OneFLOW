#include "NsCom.h"
#include "NsIdx.h"
#include "DataBase.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Chemical.h"

BeginNameSpace( ONEFLOW )

NsCom nscom;

NsCom::NsCom()
{
    init_flag = false;
}

NsCom::~NsCom()
{
    ;
}

void NsCom::Reset()
{
    init_flag = false;
    nProbe = 0;
    bcdtkey = 0;
    bcflow = nullptr;
}

void NsCom::Init()
{
    if ( init_flag ) return;
    init_flag = true;
    icmpInv = 1;
    ischeme = GetDataValue< int >( "ischeme" );
    ivischeme = GetDataValue< int >( "ivischeme" );
    timestepModel = GetDataValue< int >( "timestepModel" );
    visSRModel = GetDataValue< int >( "visSRModel" );
    visTimeStepModel = GetDataValue< int >( "visTimeStepModel" );
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
    gasModelFile = GetDataValue< std::string >( "gasModelFile" );
    machStrategy = GetDataValue< int >( "machStrategy" );

    schmidtl = GetDataValue< Real >( "schmidtl" );
    schmidtt = GetDataValue< Real >( "schmidtt" );

    max_time_ratio = GetDataValue< Real >( "max_time_ratio" );

    nEqu = GetDataValue< int >( "nEqu" );
    nTEqu = GetDataValue< int >( "nTEqu" );
    chem.Init();

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

void Extract( RealField & prim, Real & rm, Real & um, Real & vm, Real & wm, Real & pm )
{
    rm = prim[ IDX::IR ];
    um = prim[ IDX::IU ];
    vm = prim[ IDX::IV ];
    wm = prim[ IDX::IW ];
    pm = prim[ IDX::IP ];
}

bool NsCheckFunction( RealField & q )
{
    if ( q[ IDX::IR ] <= 0.0 || q[ IDX::IP ] <= 0.0 ) return false;
    return true;
}

EndNameSpace
