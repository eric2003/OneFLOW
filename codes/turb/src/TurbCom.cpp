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
#include "TurbCom.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "NsCtrl.h"
#include "DataBase.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

TurbCom turbcom;

TurbCom::TurbCom()
{
    init_flag = false;
}

TurbCom::~TurbCom()
{
    ;
}

void TurbCom::Init()
{
    if ( init_flag ) return;
    init_flag = true;
    this->nEqu = 1;
    sst_type = 0;
    if ( vis_model.vismodel == 3 )
    {
        this->nEqu = 1;
    }
    else if ( vis_model.vismodel == 4 )
    {
        this->nEqu = 2;
        if ( vis_model.visname.substr( 0, 13 ) == "2eq-kw-menter" ||
             vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            sst_type = 1;
        }
    }
    turb_ilim = GetDataValue< int >( "turb_ilim" );
    tns_ilim = GetDataValue< int >( "tns_ilim" );
    iturb_visflux = GetDataValue< int >( "iturb_visflux" );
    iprod_sa  = GetDataValue< int >( "iprod_sa" );
    iprod_sst = GetDataValue< int >( "iprod_sst" );
    easm_model = GetDataValue< int >( "easm_model" );
    ref_sa = GetDataValue< Real >( "ref_sa" );
    ref_sst = GetDataValue< Real >( "ref_sst" );
    inflow_intensity = GetDataValue< Real >( "inflow_intensity" );
    inflow_viscosity = GetDataValue< Real >( "inflow_viscosity" );
    max_vis_ratio = GetDataValue< Real >( "max_vis_ratio" );
    turb_cfl_ratio = GetDataValue< Real >( "turb_cfl_ratio" );
    rprod = GetDataValue< Real >( "rprod" );

    reynolds = GetDataValue< Real >( "reynolds" );
    oreynolds = 1.0 / reynolds;
    crelax = 1.0;

    comVis.resize( nEqu );
    flux.resize( nEqu );

    faceOuterNormal = 1;

    prims1.resize( nEqu );
    prims2.resize( nEqu );

    primt1.resize( nEqu );
    primt2.resize( nEqu );

    int nEqu_ns  = GetDataValue< int >( "nEqu" );

    ns_prims1.resize( nEqu_ns );
    ns_prims2.resize( nEqu_ns );

    q.resize( nEqu );
    q0.resize( nEqu );
    dq.resize( nEqu );

    this->inflow.resize( nEqu );

    this->InitConst();
    this->InitInflow();
}

void TurbCom::InitConst()
{
    this->pklim = 20.0;

    if ( nscom.mach_ref > 1.0 )
    {
        this->pklim = 5.0;
    }
    else
    {
        if ( iprod_sst == 0 )
        {
            // Boussinesq approximation, full production
            this->pklim = 20.0;
        }
        else if ( iprod_sst == 1 ) //SST-2003
        {
            this->pklim = 10.0;
        }
    }

    kelim = 1.0e-16;
    kwlim = 1.0e-8;
    
    a1    = 0.31;

    fbeta  = 1.0;
    fbetas = 1.0;
    betas  = 0.09;
    //实际上cmu应该为0.09，但是cmu一般都被吸收了，不需要
    cmu    = 1.0;
    clim   = 1.0;

    sigk   = 0.5;
    sigw   = 0.5;
    alphaw = 0.556;
    sigd   = 0.0;
            
    beta   = 0.83;

    sigk1   = 1.1;
    sigw1   = 0.53;
    alphaw1 = 0.518;
    beta1   = 0.0747;
    sigd1   = 1.0;
            
    sigk2   = 1.1;
    sigw2   = 1.0;
    alphaw2 = 0.44;
    beta2   = 0.0828;
    sigd2   = 0.4;

    sigma   = 2.0 / 3.0;
    osigma  = one / sigma;

    if ( vis_model.visname.substr( 0, 6 ) == "1eq-sa" )
    {
        karm   = 0.41;
        karm2  = SQR( karm );
        okarm2 = one / karm2;
        cb1    = 0.1355;
        cb2    = 0.622;
        cv1    = 7.1;
        cv13   = POWER3( cv1 );
        cw2    = 0.3;
        cw3    = 2.0;
        cw36   = pow( cw3, 6 );
        cdes   = - 1.0;
        ct3    = 1.2;
        ct4    = 0.5;
        cw1  = cb1 * okarm2 + ( one + cb2 ) / sigma;
    }
    else if ( vis_model.visname.substr( 0, 6 ) == "2eq-kw" )
    {
        if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1988" )
        {
            sigk   = 0.5;
            sigw   = 0.5;
            alphaw = 0.556;
            sigd   = 0.0;
            
            beta   = 0.83;
            betas  = 1.0;
            cmu    = 0.09;
            
            if ( easm_model == 1 )
            {
                //sigd   = 0.0 * 0.09;
                //betas  = 1.0;
                //cmu    = 0.09;
                
                //sigk   = 1.0 / 1.4;
                //sigw   = 1.0 / 2.2;
                //alphaw = 0.5615;
                //beta   = 0.83;
                //sigd   = 0.0 * 0.09;
                //betas  = 1.0;
                //cmu    = 0.09;
            }
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1993" )
        {
            sigk   = 1.0;
            sigw   = 0.6;
            alphaw = 0.556;
            beta   = 0.075;
            sigd   = 0.3;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1998" )
        {
            sigk   = 0.5;
            sigw   = 0.5;
            alphaw = 0.520;
            beta   = 0.072;
            sigd   = 0.0;
            
            if ( easm_model == 1 )
            {
                sigk   = 0.5;
                sigw   = 0.5;
                alphaw = 0.520;
                beta   = 0.83;
                sigd   = 0.0 * 0.09;
                betas  = 1.0;
                cmu    = 0.09;
            }
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-2006" )
        {
            sigk   = 0.6;
            sigw   = 0.5;
            alphaw = 0.520;
            beta   = 0.0708;
            sigd   = 0.0;
            clim   = 0.875;
        }
        else if ( vis_model.visname.substr( 0, 14 ) == "2eq-kw-kok-tnt" )
        {
            //sigk   = 0.5;
            sigk   = 0.667;
            sigw   = 0.5;
            alphaw = 0.556;
            beta   = 0.075;
            sigd   = 0.5;
            
            if ( easm_model == 1 )
            {
                sigk   = 1.0 / 1.4;
                sigw   = 1.0 / 2.2;
                alphaw = 0.5615;
                beta   = 0.83;
                sigd   = 0.5 * 0.09;
                betas  = 1.0;
                cmu    = 0.09;
            }
        }
        else if ( vis_model.visname.substr( 0, 17 ) == "2eq-kw-menter-sst" )
        {
            sigk1   = 0.85;
            sigw1   = 0.5;
            alphaw1 = 0.5532;
            beta1   = 0.075;
            
            sigk2   = 1.0;
            sigw2   = 0.856;
            alphaw2 = 0.4404;
            beta2   = 0.0828;
        }
        else if ( vis_model.visname.substr( 0, 17 ) == "2eq-kw-menter-bsl" )
        {
            sigk1   = 0.5;
            sigw1   = 0.5;
            alphaw1 = 0.5532;
            beta1   = 0.075;
            
            sigk2   = 1.0;
            sigw2   = 0.856;
            alphaw2 = 0.4404;
            beta2   = 0.0828;
        }
        else if ( vis_model.visname.substr( 0, 10 ) == "2eq-kw-hyb" )
        {
        }
        else if ( vis_model.visname.substr( 0, 10 ) == "2eq-kw-bdp" )
        {
        }
    }
    else if ( vis_model.visname.substr( 0, 4 ) == "easm" )
    {
        // 专门的easm，不是随意组合的那种
        if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2001" )
        {
            sigk   = 0.5;
            sigw   = 0.5339;
            alphaw = 0.575;
            beta   = 0.83;
            sigd   = 0.0;
            betas  = 1.0;
            cmu    = 0.0895;
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2003" )
        {
            sigk   = 1.0;
            sigw   = 0.5339;
            alphaw = 0.53;
            beta   = 0.83;
            sigd   = 0.0;
            betas  = 1.0;
            cmu    = 0.0895;
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            sigk1   = 1.1;
            sigw1   = 0.53;
            alphaw1 = 0.518;
            beta1   = 0.0747;
            sigd1   = 1.0;
            
            sigk2   = 1.1;
            sigw2   = 1.0;
            alphaw2 = 0.44;
            beta2   = 0.0828;
            sigd2   = 0.4;
            
            betas  = 0.09;
        }
    }


    sac2 = 0.7;
    sac3 = 0.9;

    cb2s = osigma * cb2;
    cw1k = cw1 * karm2;

    fwStar  = 0.424;
    cdes    = 0.65;
    cdes_ke = 0.61;
    cdes_kw = 0.78;
}

void TurbCom::InitInflow()
{
    Real ref_density = nscom.inflow[ IDX::IR ];
    if ( this->nEqu == 1 )
    {
        if ( vis_model.vismodel == 3 )
        {
            inflow[ ISA ] = ref_sa;
        }
        else
        {
            inflow[ ISA ] = 0.1 / ( ref_density * cmu + SMALL );
        }
    }
    else if ( this->nEqu >= 2 )
    {
        Real kwoo = 1.0;
        Real mutoo = 0.001;

        if ( vis_model.vismodel == 4 )
        {
            mutoo = ref_sst;
        }

        Real keoo = mutoo * kwoo / ( cmu * reynolds * ref_density );
        inflow[ IKE ] = keoo;
        inflow[ IKW ] = kwoo;
    }

    if ( transition_model == ITReGama ) 
    {
        trans.ce1 = 1.0;
        trans.ca1 = 2.0;
        trans.ce2 = 50.0;
        trans.ca2 = 0.06;
        trans.dct = 2.0;
        trans.df  = 1.0;
        trans.cct = 0.03;
        trans.s1  = 2.0;

        Real kwoo = 1.0;
        Real mutoo = inflow_viscosity;
        Real keoo = mutoo * kwoo / ( cmu * reynolds * ref_density );
        
        inflow[ IKE ] = keoo;
        inflow[ IKW ] = kwoo;
        
        inflow[ ITGama ] = 1.0;

        Real TUoo = 0.18;
        Real Fcta = 1.0;

        if ( inflow_intensity <= 0.0 )
        {
            TUoo = trans.CalcIntensity( 1.0, keoo );
        }
        else
        {
            keoo = 1.5 * SQR( 0.01 * inflow_intensity );
            kwoo = keoo * ( cmu * reynolds * ref_density ) / mutoo;
            inflow[ IKE ] = keoo;
            inflow[ IKW ] = kwoo;
        }
        
        inflow[ ITRect ] = trans.EmpiricalCorrelationOfRectat( TUoo, Fcta );
    }

}

void TurbCom::CalcSigkw()
{
    if ( sst_type == 0 ) return;

    bld = half * ( bld1 + bld2 );
    sigk = bld * sigk1 + ( 1.0 - bld ) * sigk2;
    sigw = bld * sigw1 + ( 1.0 - bld ) * sigw2;
}

void TurbCom::CalcWorkVar()
{
    //work1 = ui,i
    //work2 = u1,1^2 + u2,2^2 + u3,3^2
    //work3 = ( u1,2+u2,1 )^2 + ( u1,3+u3,1 )^2 + ( u2,3+u3,2 )^2
    //work4 = sqrt( (u3,2-u2,3)^2 + (u1,3-u3,1)^2+ ( u2,1-u1,2 )^2 )
    work1 = dudx + dvdy + dwdz;
    work2 = SQR( dudx, dudy, dudz );
    work3 = SQR( dudy + dvdx,  dvdz + dwdy, dwdx + dudz );
    work4 = DIST( dwdy - dvdz, dudz - dwdx, dvdx - dudy );
}

void TurbCom::CalcSaProd()
{
    if ( iprod_sa == 1 )
    {
        //classical production ( instead of rotational )
        //stp = ( 2 * ( ui,i^2 ) + ( ui,j + uj,i )^2 - 2/3 * ( ui,i )^2 )
        stp   = two * work2 + work3 - two3rd * work1 * work1;
        work4 = sqrt( ABS( stp ) );
    }
    else if ( iprod_sa == 4 )
    {
        //correction on the rotational production term for laminar vortex
        stp   = two * work2 + work3 - two3rd * work1 * work1;
        str   = work4;
        work4 = str - two * MAX( zero, str - sqrt( ABS( stp ) ) );
    }
    else if ( iprod_sa == 5 )
    {
        //higher correction on the rotational production term for laminar vortex
        stp   = two * work2 + work3 - two3rd * work1 * work1;
        str   = work4;
        stp   = ABS( stp );
        work4 = str * ( one - MAX( zero , ( str * str - stp ) / ( str * str + stp + SMALL ) ) );
    }
}

void TurbCom::CalcVGrad()
{
    s11 = dudx;
    s22 = dvdy;
    s33 = dwdz;
    s12 = half * ( dudy + dvdx );
    s13 = half * ( dudz + dwdx );
    s23 = half * ( dvdz + dwdy );
    w12 = half * ( dudy - dvdx );
    w13 = half * ( dudz - dwdx );
    w23 = half * ( dvdz - dwdy );
                                           
    sij2 = two * ( SQR( s11, s22, s33 ) + two * SQR( s12, s13, s23 ) );
    divv = s11 + s22 + s33;
}

void TurbCom::CalcProdk()
{
    if ( iprod_sst == 0 )
    {
        // Boussinesq approximation, full production
        prodk = vist * ( sij2 - two3rd * SQR( divv ) ) * oreynolds - two3rd * rho * ke * divv;
    }
    else if ( iprod_sst == 1 ) //SST-2003
    {
        prodk = vist * sij2 * oreynolds;
    }
    else //Vorticity Source Term (SST-V) 
    {
        Real vort2 = four * SQR( w12, w13, w23 );
        prodk = vist * ( vort2 ) * oreynolds - two3rd * rho * ke * divv;
    }
}

void TurbCom::CalcDissk()
{
    if ( des_model ) 
    {
        dissk = fbetas * rho * ke * sqrt( ke ) / len_scale;
    }
    else
    {
        dissk = fbetas * betas * rho * ke * kw;
    } 
}

void TurbCom::LimitProdk()
{
    prodk = MIN( prodk, pklim * dissk );
}

void TurbCom::CalcProdwKwMenter()
{
    beta   = bld * beta1  + ( 1.0 - bld ) * beta2;
    alphaw = bld * alphaw1 + ( 1.0 - bld ) * alphaw2;

    prodw = alphaw * rho / ( vist + SMALL ) * prodk * reynolds;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = two * ( one - bld ) * rho * sigw2 * cross_term / ( kw + SMALL );
}

void TurbCom::CalcProdwKwWilcox1998()
{
    prodw = alphaw * kw / ( ke + SMALL ) * prodk;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = sigd * rho * MAX( cross_term, zero ) / ( kw + SMALL );
}

void TurbCom::CalcProdwKwWilcox2006()
{
    sigd = zero;
    if ( cross_term > zero )
    {
        // 1/8
        sigd = 0.125;
    }

    prodw = alphaw * kw / ( ke + SMALL ) * prodk;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = sigd * rho * MAX( cross_term, zero ) / ( kw + SMALL );
}

void TurbCom::CalcProdwKwDefault()
{
    prodw = alphaw * kw / ( ke + SMALL ) * prodk;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = sigd * rho * MAX( cross_term, zero ) / ( kw + SMALL );
}

void TurbCom::CalcProdwEasmKw2003()
{
    prodw = alphaw * kw / ( ke + SMALL ) * prodk;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = sigd * rho * MAX( cross_term, zero ) / ( kw + SMALL );
}

void TurbCom::CalcProdwEasmKw2005()
{
    beta   = bld * beta1   + ( 1.0 - bld ) * beta2;
    alphaw = bld * alphaw1 + ( 1.0 - bld ) * alphaw2;
    sigd   = bld * sigd1   + ( 1.0 - bld ) * sigd2;

    prodw = alphaw * kw / ( ke + SMALL ) * prodk;
    dissw = fbeta * beta * rho * SQR( kw );
    cdkww = sigd * rho * MAX( cross_term, zero ) / ( kw + SMALL );
}

void TurbCom::ModifyPd()
{
    if ( transition_model == ITReGama ) 
    {                    
        prodk = turb_trans.CorrectionOfProductionInKEquation ( gmeff, prodk );
        dissk = turb_trans.CorrectionOfDestructionInKEquation( gmeff, dissk );
    }
}

void TurbCom::CalcSrc()
{
    srck  = prodk - dissk;
    srcw  = prodw - dissw + cdkww;

    // The linearization proposed by Menter is used.
    diak = - two * fbetas * betas * kw;
    diaw = - two * fbeta  * beta  * kw - ABS( cdkww ) / ( rho * kw + SMALL );

    fskn = half * ( srck - ABS( srck ) );
    fswn = half * ( srcw - ABS( srcw ) );

    oork = 1.0 / ( rho * ke + SMALL );
    oorw = 1.0 / ( rho * kw + SMALL );
}

void TurbCom::CalcCrossDiff()
{
    crossdiff = 2.0 * rho * cross_term * sigw2 / ( kw + SMALL );
    if ( crossdiff > cdkwmax )
    {
        cdkwmax = crossdiff;
    }
}

void TurbCom::CalcCdkwmin()
{
    //cdkwmin = cdkwmax * 1.0e-8;
    cdkwmin = 1.0e-20;
}

void TurbCom::CalcCellBlendingTerm()
{
    //calculate cd_kw
    crossdiff = 2.0 * rho * sigw2 * cross_term / ( kw + SMALL );
    Real cdkw;

    if ( iprod_sst == 1 )
    {
        cdkw = MAX( crossdiff, 1.0e-10 );
    }
    else //Original Menter CD_kw calculation
    {
        cdkw = MAX( crossdiff, 1.0e-20 );
    }

    Real dist2  = SQR( dist );

    //calculate arg1
    Real term1 = sqrt( ABS( ke ) )/( betas * dist * kw + SMALL );
    Real term2 = 500.0 * visl / ( rho * kw * dist2 * reynolds + SMALL );
    Real term3 = MAX( term1, term2 );
    Real term4 = 4.0 * rho * sigw2 * ke / ( cdkw * dist2 + SMALL );
    Real arg1  = MIN( term3, term4 );

    bld = tanh( POWER4( arg1 ) );
}

void TurbCom::CalcCrossing()
{
    cross_term = dkedx * dkwdx + dkedy * dkwdy + dkedz * dkwdz;
}

void TurbCom::CalcFbetaOfKwWilcox1998()
{
    xk = cross_term / ( pow( kw, 3 ) + SMALL );
    if ( xk <= 0.0 )
    {
        fk = 1.0;
    }
    else
    {
        xk2 = SQR( xk );
        // Bounded crossing_term diffusion function( 1<fk<1.7 )
        fk  = ( one + 680.0 * xk2 ) / ( one + 400.0 * xk2 );
    }
                
    xw = - ( w12 * w12 * ( s11 + s22 ) + w13 * w13 * ( s11 + s33 ) + w23 * w23 * ( s22 + s33 ) )
         + two * ( - w13 * w23 * s12 + w12 * w23 * s13 - w12 * w13 * s23 );
    xw = ABS( xw /( pow( betas * kw, 3 ) + SMALL ) );
    // Bounded vortex stretching function( 7/8<fw<1 )
    fw = ( one + 70.0 * xw ) / ( one + 80.0 * xw );
                
    fbetas = fk;
    fbeta  = fw;
}

void TurbCom::CalcFbetaOfKwWilcox2006()
{
    Real sh11 = s11 - half * divv;
    Real sh22 = s22 - half * divv;
    Real sh33 = s33 - half * divv;
    Real sh12 = s12;
    Real sh13 = s13;
    Real sh23 = s23;
                
    fk = 1.0;
    xw = - ( w12 * w12 * ( sh11 + sh22 ) + w13 * w13 * ( sh11 + sh33 ) + w23 * w23 * ( sh22 + sh33 ) )
         + two * ( - w13 * w23 * sh12 + w12 * w23 * sh13 - w12 * w13 * sh23 );
    xw = ABS( xw /( pow( betas * kw, 3 ) + SMALL ) );
    // Bounded vortex stretching function( 7/8<fw<1 )
    fw = ( one + 85 * xw ) / ( one + 100.0 * xw );
                
    fbetas = fk;
    fbeta  = fw;
}

void TurbCom::CalcFbetaOfEasmKw2003()
{
    xk = cmu * cmu * cross_term / ( pow( kw, 3 ) + SMALL );
    fk = 1.0;
    if ( xk > 0.0 )
    {
        xk2 = xk * xk;
        // Bounded crossingTermField diffusion function( 1<fk<1.7 )
        fk  = ( one + 680.0 * xk2 )/( one + 400.0 * xk2 );
    }

    fw = 1.0;
                        
    fbetas = fk;
    fbeta  = fw;
}

void TurbCom::RGamaTransition()
{
    if ( transition_model != ITReGama ) return;

    Real vorx   = dvdz - dwdy;
    Real vory   = dudz - dwdx;
    Real vorz   = dudy - dvdx;
                            
    Real vorticity  = DIST( vorx, vory, vorz );
    Real strainRate = sqrt( sij2 );
                            
    Real absU    = MAX( DIST( um, vm, wm ), SMALL );
                            
    Real RT      = turb_trans.ViscosityRatio( rho, visl, ke, kw, reynolds );
    Real Rev     = turb_trans.ReynoldsNumberBasedOnStrainRate( rho, dist, visl, strainRate, reynolds );
    Real Rectac  = turb_trans.TransitionOnsetMomentumThicknessReynolds( rectabar );

    Real Rew     = turb_trans.ReynoldsNumberBasedOnDissipation( rho, rectabar, visl, kw, reynolds );
    Real Flength = turb_trans.HighReynoldsCorrectionOfFlength( Rew, turb_trans.FlengthGivenByLangtry( rectabar ) );
    Real Fonset  = turb_trans.TransitionOnsetFunction( Rev, Rectac, RT );
    Real Fturb   = turb_trans.ControlFunctionFturb( RT );
                            
    Real production1OfGama = turb_trans.ca1 * Flength * rho * strainRate * sqrt( intermittency * Fonset );
    Real production2OfGama = production1OfGama * ( - turb_trans.ce1 * intermittency );
                            
    Real productionOfGama  = production1OfGama + production2OfGama;
                            
    Real specGamap = MIN( 0.0, 1.5 * production2OfGama + 0.5 * production1OfGama ) / ( rho * intermittency );
                            
    Real destruction2OfGama = - turb_trans.ca2 * Fturb * rho * vorticity * intermittency;
    Real destruction1OfGama = - turb_trans.ce2 * intermittency * destruction2OfGama;
                            
    Real destructionOfGama = destruction1OfGama + destruction2OfGama;
    Real specGamad = MAX( 0.0, 2.0 * destruction1OfGama + destruction2OfGama ) / ( rho * intermittency );

    srcg = ( productionOfGama - destructionOfGama );
    speg = ( specGamad - specGamap );
                            
    Real Fctat = turb_trans.BlendingFunctionOfFctat( intermittency, Rew, rectabar, vorticity, visl, rho, absU, dist, reynolds );
    Real tscl  = turb_trans.TimeScaleInSourceTerm(  rho, absU, visl, reynolds ); 
    Real TU    = turb_trans.CalcIntensity( absU, ke );
    gmeff = turb_trans.SeparationCorrectionOfIntermittency( intermittency, Rev, Rectac, RT, Fctat );
                            
    Real dUds  = turb_trans.AccelerationAlongStreamline( um, vm, wm, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz );
                            
    Real momentumThickness = 0.0;
    Real Rectat;
    for ( int iter = 0; iter < 10; ++ iter )
    {
        Real lamdacta  = turb_trans.PressureGradientFunction( rho, momentumThickness, visl, dUds, reynolds );
        Real Flamdacta = turb_trans.EmpiricalCorrelationOfFlamdacta( TU, lamdacta );
        Rectat = turb_trans.EmpiricalCorrelationOfRectat( TU, Flamdacta );
        momentumThickness = turb_trans.MomentumThickness( rho, absU, visl, Rectat, reynolds );
    }
                            
    Real coecommon = turb_trans.cct * ( 1.0 - Fctat ) / MAX( tscl, 1.0e-20 );
                            
    Real prodRe    = coecommon * rho * Rectat;
    Real destRe    = coecommon * rho * rectabar;

    srcr = ( prodRe - destRe );
    sper = destRe / ( rho * rectabar );
}

void TurbCom::CalcSrcSa()
{
    Real olam = rho / ( visl + SMALL );
    Real d2   = SQR( len_scale );
    Real od2  = one / d2;

    CalcWorkVar();
    CalcSaProd();
    this->str = this->work4;

    //this->compress = this->nuet * this->work1;
    this->compress = 0.0;

    xsi  = nuet * olam + SMALL;
    xsi3 = POWER3( xsi );
    fv1  = xsi3 / ( xsi3 + cv13 );
    fv2  = one - xsi / ( one + xsi * fv1 );

    //here we absorb oreynolds into rs
    Real rs = nuet * okarm2 * od2 * oreynolds;
    Real nuetRs = nuet * rs;

    Real sBar    = rs * fv2;
    Real omega   = str;

    Real std;

    if ( sBar >= - sac2 * omega )
    {
        std = omega + sBar;
    }
    else
    {
        Real term1 = ( sac2 * sac2 * omega + sac3 * sBar );
        Real term2 = ( sac3 - 2 * sac2 ) * omega - sBar;
        std = omega + omega * term1 / term2;
    }

    Real ostd   = one / ( std + SMALL );
    Real r      = rs * ostd;

    r    = MIN( r, ten );
    Real r5   = pow( r, 5 );
    Real g    = r + cw2 * ( r5 * r - r );
    Real g6   = pow( g, 6 );
    Real or6  = one / six;
    fw   = g * pow( ( one + cw36 ) / ( g6 + cw36 ), or6 );

    //calculate "negative" production term in order to assure a zero
    //solution in the laminar region
    ftrans  = zero;
    xsi2 = SQR( xsi );
    ft2  = 0.0;

    if ( ft2_flag )
    {
        ftrans = one;
        ft2 = ct3 * exp( - ct4 * xsi2 ) * ftrans;
    }

    Real grd2 = SQR( dqdxSa, dqdySa, dqdzSa );

    Real prod = cb1 * ( one - ft2 ) * std * nuet;    //cb1  * omega * muet + cb1 * muetRs * fv2;
    Real diff = cb2s * grd2 * oreynolds;             //cb2s * density * gradnue2;
    Real dest = ( cw1k * fw - cb1 * ft2 ) * nuetRs;     //cw1k * muetRs * fw

    diff = MIN( rprod * prod, diff );

    Real srcTerm  = prod - dest;

    srcTerm += diff; 
        
    Real dfv2dk = ( three * cv13 * pow( fv1 / ( xsi + SMALL ), 2 ) - one ) / pow( one + xsi * fv1, 2 );
    Real dsdnu  = oreynolds * okarm2 * od2 * ( fv2 + xsi * dfv2dk );

    Real dfwdg   = fw / ( g + SMALL ) * ( one - g6 / ( g6 + cw36 ) );
    Real dgdr    = one + cw2 * ( six * r5 - one );
    Real drdnu   = oreynolds * okarm2 * od2 * ostd * ( one - nuet * ostd * dsdnu );
    Real dfwdnu  = dfwdg * dgdr * drdnu;

    Real dft2dnu = 0.0;

    if ( ABS( ft2 ) > 1.0e-5 )
    {
        dft2dnu = - two * ct3 * ct4 * xsi * exp( - ct4 * xsi2 ) * olam;
    }

    Real prod0 = cb1 * ( one - ft2 ) * std;
    Real dest0 = ( cw1 * fw - cb1 * ft2 * okarm2 ) * nuet * od2 * oreynolds;

    Real prodp = cb1 * ( ( one - ft2 ) * dsdnu );
    Real destp = ( cw1 * ( fw + dfwdnu * nuet ) - cb1 * ft2 * okarm2 ) * od2 * oreynolds;

    Real prde0 = prod0 - dest0;
    Real prdep = prodp - destp;

    Real part1 = - half * ( prde0 - ABS( prde0 ) );
    Real part2 = - half * ( prdep - ABS( prdep ) ) * nuet;

    spec_sa = part1 + part2;
    res_sa  = srcTerm;
}

EndNameSpace