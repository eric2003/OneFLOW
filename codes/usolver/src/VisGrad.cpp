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

#include "VisGrad.h"
#include "Com.h"
#include "HXMath.h"
#include "Boundary.h"
#include "UCom.h"
#include "BcRecord.h"

BeginNameSpace( ONEFLOW )

VisGradGeom vgg;

VisGradGeom::VisGradGeom()
{
    skewAngle = 20.0;
}

VisGradGeom::~VisGradGeom()
{
    ;
}

void VisGradGeom::CalcFaceWeight()
{
    dxl = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
    dyl = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
    dzl = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

    dxr = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
    dyr = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
    dzr = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

    delt1 = DIST( dxl, dyl, dzl );
    delt2 = DIST( dxr, dyr, dzr );
    delta = 1.0 / ( delt1 + delt2 + SMALL );

    fw1 = delt2 * delta;
    fw2 = delt1 * delta;
}

void VisGradGeom::CalcAngle( Real dx, Real dy, Real dz, Real dist, Real & angle )
{
    Real result = dist / ( DIST( dx, dy, dz ) + SMALL );
    if ( result >   1.0 ) result =   1.0;
    if ( result < - 1.0 ) result = - 1.0;
    angle = asin( result ) * 180.0 / PI;
}

void VisGradGeom::PrepareCellGeom()
{
    this->dxl = ( * ug.xcc )[ ug.lc ] - ( * ug.xfc )[ ug.fId ];
    this->dyl = ( * ug.ycc )[ ug.lc ] - ( * ug.yfc )[ ug.fId ];
    this->dzl = ( * ug.zcc )[ ug.lc ] - ( * ug.zfc )[ ug.fId ];

    this->dxr = ( * ug.xcc )[ ug.rc ] - ( * ug.xfc )[ ug.fId ];
    this->dyr = ( * ug.ycc )[ ug.rc ] - ( * ug.yfc )[ ug.fId ];
    this->dzr = ( * ug.zcc )[ ug.rc ] - ( * ug.zfc )[ ug.fId ];

    this->d1  = gcom.xfn * this->dxl + gcom.yfn * this->dyl + gcom.zfn * this->dzl;
    this->d2  = gcom.xfn * this->dxr + gcom.yfn * this->dyr + gcom.zfn * this->dzr;

    this->dxnl  = gcom.xfn * this->d1 - this->dxl;
    this->dynl  = gcom.yfn * this->d1 - this->dyl;
    this->dznl  = gcom.zfn * this->d1 - this->dzl;

    this->dxnr  = gcom.xfn * this->d2 - this->dxr;
    this->dynr  = gcom.yfn * this->d2 - this->dyr;
    this->dznr  = gcom.zfn * this->d2 - this->dzr;

    //d1=( d1x, d1y ) = d1 * ( fnx, fny )
    //( r1x, r1y ) - ( rLx, rLy ) = ( dr1x, dr1y ) - ( drLx, drLy )
    //( dr1x, dr1y ) - ( drLx, drLy ) = d1 * ( fnx, fny ) - ( dxl, dyl )

    this->CalcAngle( this->dxl, this->dyl, this->dzl, - this->d1, this->angle1 );
    this->CalcAngle( this->dxr, this->dyr, this->dzr,   this->d2, this->angle2 );
}

void VisGradGeom::CalcGradCoef()
{
    this->dx  = ( * ug.xcc )[ ug.rc ] - ( * ug.xcc )[ ug.lc ];
    this->dy  = ( * ug.ycc )[ ug.rc ] - ( * ug.ycc )[ ug.lc ];
    this->dz  = ( * ug.zcc )[ ug.rc ] - ( * ug.zcc )[ ug.lc ];

    this->ods = 1.0 / DIST( this->dx, this->dy, this->dz );

    this->dx *= this->ods;
    this->dy *= this->ods;
    this->dz *= this->ods;
}


VisGrad::VisGrad()
{
    ;
}

VisGrad::~VisGrad()
{
    ;
}

void VisGrad::Init( int nEqu )
{
    this->nEqu = nEqu;
    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    q11.resize( nEqu );
    q22.resize( nEqu );

    dqdx.resize( nEqu );
    dqdy.resize( nEqu );
    dqdz.resize( nEqu );
    dqdn.resize( nEqu );

    dqdx1.resize( nEqu );
    dqdy1.resize( nEqu );
    dqdz1.resize( nEqu );

    dqdx2.resize( nEqu );
    dqdy2.resize( nEqu );
    dqdz2.resize( nEqu );

    dqdn1.resize( nEqu );
    dqdn2.resize( nEqu );

    dqdt1.resize( nEqu );
    dqdt2.resize( nEqu );
}

void VisGrad::AverGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdx[ iEqu ] = half * ( dqdx1[ iEqu ] + dqdx2[ iEqu ] );
        dqdy[ iEqu ] = half * ( dqdy1[ iEqu ] + dqdy2[ iEqu ] );
        dqdz[ iEqu ] = half * ( dqdz1[ iEqu ] + dqdz2[ iEqu ] );
    }
}

void VisGrad::ZeroNormalGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn[ iEqu ] = 0.0;
    }
}

void VisGrad::AverFaceValue()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        q[ iEqu ] = half * ( q1[ iEqu ] + q2[ iEqu ] );
    }
}

void VisGrad::CorrectFaceGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        CorrectGrad( q1[ iEqu ], q2[ iEqu ], dqdx[ iEqu ], dqdy[ iEqu ], dqdz[ iEqu ], vgg.dx, vgg.dy, vgg.dz, vgg.ods );
    }
}

void VisGrad::CalcNormalGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn[ iEqu ] = gcom.xfn * dqdx[ iEqu ] + gcom.yfn * dqdy[ iEqu ] + gcom.zfn * dqdz[ iEqu ];
    }
}

bool VisGrad::FaceAngleIsValid()
{
    // Theoretically, more accurate to include the following terms
    bool result =  vgg.angle1 > vgg.skewAngle && vgg.angle2 > vgg.skewAngle;
    return result;
}

bool VisGrad::TestSatisfied()
{
    bool result = vgg.angle1 > 0.0 && vgg.angle2 > 0.0 && ABS( vgg.d1 ) > SMALL && ABS( vgg.d2 ) > SMALL;
    if ( result )
    {
        this->CalcC1C2();
    }
    return result;
}

bool VisGrad::New1Satisfied()
{
    bool result =  vgg.d1 * vgg.d2 < 0.0 && ABS( vgg.d1 ) > SMALL && ABS( vgg.d2 ) > SMALL;
    if ( result )
    {
        vgg.d  = - two * vgg.d1 * vgg.d2 / ( SQR( vgg.d1, vgg.d2 ) + SMALL );
        vgg.od = vgg.d / ( vgg.d2 - vgg.d1 );
    }

    return result;
}

bool VisGrad::New2Satisfied()
{
    vgg.d = - two *  vgg.d1 *  vgg.d2 / ( SQR(  vgg.d1,  vgg.d2 ) + SMALL );

    bool result =  vgg.d > 0.1;

    return result;
}

void VisGrad::CalcC1C2()
{
    Real dtmp = SQR( vgg.d1, vgg.d2 );
    vgg.c1 = SQR( vgg.d1 ) / dtmp;
    vgg.c2 = 1.0 - vgg.c1;

    if ( ug.fId < ug.nBFace )
    {
        int bcType = ug.bcRecord->bcType[ ug.fId ];

        if ( bcType != BC::INTERFACE )
        {
            vgg.c1 = 1.0;
            vgg.c2 = 0.0;
        }
    }
}

void VisGrad::AccurateSideValue()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        q11[ iEqu ] += dqdx1[ iEqu ] * vgg.dxnl + dqdy1[ iEqu ] * vgg.dynl + dqdz1[ iEqu ] * vgg.dznl;
        q22[ iEqu ] += dqdx2[ iEqu ] * vgg.dxnr + dqdy2[ iEqu ] * vgg.dynr + dqdz2[ iEqu ] * vgg.dznr;
    }
}

void VisGrad::AccurateFaceValue()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        q[ iEqu ] = half * ( q1[ iEqu ] + q2[ iEqu ] );
    }
}

void VisGrad::ModifyFaceGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdx[ iEqu ] = vgg.fw1 * dqdx1[ iEqu ] + vgg.fw2 * dqdx2[ iEqu ];
        dqdy[ iEqu ] = vgg.fw1 * dqdy1[ iEqu ] + vgg.fw2 * dqdy2[ iEqu ];
        dqdz[ iEqu ] = vgg.fw1 * dqdz1[ iEqu ] + vgg.fw2 * dqdz2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdt1[ iEqu ] = gcom.t1x * dqdx[ iEqu ] + gcom.t1y * dqdy[ iEqu ] + gcom.t1z * dqdz[ iEqu ];
        dqdt2[ iEqu ] = gcom.t2x * dqdx[ iEqu ] + gcom.t2y * dqdy[ iEqu ] + gcom.t2z * dqdz[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdx[ iEqu ] = gcom.xfn * dqdn[ iEqu ] + gcom.t1x * dqdt1[ iEqu ] + gcom.t2x * dqdt2[ iEqu ];
        dqdy[ iEqu ] = gcom.yfn * dqdn[ iEqu ] + gcom.t1y * dqdt1[ iEqu ] + gcom.t2y * dqdt2[ iEqu ];
        dqdz[ iEqu ] = gcom.zfn * dqdn[ iEqu ] + gcom.t1z * dqdt1[ iEqu ] + gcom.t2z * dqdt2[ iEqu ];
    }
}

void VisGrad::CalcTestMethod()
{
    if ( this->FaceAngleIsValid() )
    {
        this->AccurateSideValue();
        this->AccurateFaceValue();
    }

    if ( ! this->TestSatisfied() ) return;

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn1[ iEqu ] = ( q11[ iEqu ] - q[ iEqu ] ) / vgg.d1;
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn2[ iEqu ] = ( q22[ iEqu ] - q[ iEqu ] ) / vgg.d2;
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn[ iEqu ] = vgg.c1 * dqdn1[ iEqu ] + vgg.c2 * dqdn2[ iEqu ];
    }
}

void VisGrad::CalcNew1Method()
{
    if ( ! this->New1Satisfied() ) return;

    this->AccurateSideValue();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn[ iEqu ] = ( q22[ iEqu ] - q11[ iEqu ] ) * vgg.od;
    }
}

void VisGrad::CalcNew2Method()
{
    if ( ! this->New2Satisfied() ) return;

    // Theoretically, more accurate to include the following terms
    this->AccurateSideValue();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn1[ iEqu ] = ( q11[ iEqu ] - q[ iEqu ] ) / vgg.d1;
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dqdn2[ iEqu ] = ( q22[ iEqu ] - q[ iEqu ] ) / vgg.d2;
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        if ( dqdn1[ iEqu ] * dqdn2[ iEqu ] > 0.0 )
        {
            dqdn[ iEqu ] = two * dqdn1[ iEqu ] * dqdn2[ iEqu ] / ( dqdn1[ iEqu ] + dqdn2[ iEqu ] );
        }
    }
}

void CorrectGrad( Real fl, Real fr, Real & dfdx, Real & dfdy, Real & dfdz, Real dx, Real dy, Real dz, Real ods )
{
    Real fOld, fnew, df;
    fOld = dfdx * dx + dfdy * dy + dfdz * dz;
    fnew = ( fr - fl ) * ods;
    df   = fnew - fOld;

    dfdx += df * dx;
    dfdy += df * dy;
    dfdz += df * dz;
}

EndNameSpace