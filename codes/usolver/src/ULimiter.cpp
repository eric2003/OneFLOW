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

#include "ULimiter.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "Ctrl.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "Iteration.h"

BeginNameSpace( ONEFLOW )

Real BarthFunction( Real dq, Real minv, Real maxv, Real dot )
{
    if ( dot < 0.0 ) 
    {
        return 0.0;
    }
    else if ( dq > maxv )
    {
        return maxv / dq;
    }
    else if ( dq < minv )
    {
        return minv / dq;
    }

    return 1.0;
}

Real Vencat( Real x, Real y, Real eps )
{
    Real x2 = x * x;
    Real xy = x * y;
    Real y2 = y * y;
    Real cc = 2.0;
    //Real cc = 0.1;
    return ( ( x2 + eps + cc * xy ) / ( x2 + cc * y2 + xy + eps ) );
}

Real VencatEric( Real x, Real y, Real eps )
{
    Real x2 = x * x;
    Real xy = x * y;
    Real y2 = y * y;
    return ( ( x2 + eps + 1.0 * xy ) / ( x2 + 2.0 * y2 + xy + eps ) );
}

Real VencatFunction( Real dq, Real minv, Real maxv, Real eps )
{
    if ( dq > SMALL )
    {
        return Vencat( maxv, dq, eps );
    }
    else if ( dq < - SMALL )
    {
        return Vencat( minv, dq, eps );
    }

    return 1.0;
}

Real VencatFunctionEric( Real dq, Real minv, Real maxv, Real eps )
{
    if ( dq > SMALL )
    {
        return VencatEric( maxv, dq, eps );
    }
    else if ( dq < - SMALL )
    {
        return VencatEric( minv, dq, eps );
    }

    return 1.0;
}

Lim::Lim()
{
    ;
}

Lim::~Lim()
{
    ;
}

LimField::LimField()
{
    ;
}

LimField::~LimField()
{
    ;
}

Real LimField::ModifyLimiter( Real phil, Real phir )
{
    if ( ug.ireconface == 0 )
    {
        return phil;
    }
    else
    {
        return MIN( phil, phir );
    }
}

void LimField::GetQlQr()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 433 )
        {
            vector< Real > tmp1, tmp2;
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                tmp1.push_back( ( * this->q )[ iEqu ][ ug.lc ] );
                tmp2.push_back( ( * this->q )[ iEqu ][ ug.rc ] );
            }
            int kkk = 1;
        }

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            ( * this->qf1 )[ iEqu ][ ug.fId ] = ( * this->q )[ iEqu ][ ug.lc ];
            ( * this->qf2 )[ iEqu ][ ug.fId ] = ( * this->q )[ iEqu ][ ug.rc ];
        }
    }
}

void LimField::BcQlQrFix()
{
    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        int bcType = ug.bcRecord->bcType[ fId ];
        if ( bcType == BC::INTERFACE ) continue;
        if ( bcType == BC::PERIODIC  ) continue;

        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real tmp = half * ( ( * this->q )[ iEqu ][ ug.lc ] + ( * this->q )[ iEqu ][ ug.rc ] );

            ( * this->qf1 )[ iEqu ][ ug.fId ] = tmp;
            ( * this->qf2 )[ iEqu ][ ug.fId ] = tmp;
        }
    }
}


void LimField::CalcFaceValue()
{
    RealField qTry( this->nEqu );

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 433 )
        {
            int kkk = 1;
        }

        Real dx = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        Real dy = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        Real dz = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            qTry[ iEqu ] = ( * this->qf1 )[ iEqu ][ ug.fId ];
        }

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real dqdx = ( * this->dqdx )[ iEqu ][ ug.lc ];
            Real dqdy = ( * this->dqdy )[ iEqu ][ ug.lc ];
            Real dqdz = ( * this->dqdz )[ iEqu ][ ug.lc ];

            Real phil  = ( * this->limiter )[ iEqu ][ ug.lc ];
            Real phir  = ( * this->limiter )[ iEqu ][ ug.rc ];
            Real phi = this->ModifyLimiter( phil, phir );

            qTry[ iEqu ] += phi * ( dqdx * dx + dqdy * dy + dqdz * dz );
        }

        if ( ( * this->ckfun )( qTry ) )
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf1 )[ iEqu ][ ug.fId ] = qTry[ iEqu ];
            }
        }

        dx = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        dy = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        dz = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            qTry[ iEqu ] = ( * this->qf2 )[ iEqu ][ ug.fId ];
        }

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real dqdx = ( * this->dqdx )[ iEqu ][ ug.rc ];
            Real dqdy = ( * this->dqdy )[ iEqu ][ ug.rc ];
            Real dqdz = ( * this->dqdz )[ iEqu ][ ug.rc ];

            Real phil = ( * this->limiter )[ iEqu ][ ug.lc ];
            Real phir = ( * this->limiter )[ iEqu ][ ug.rc ];
            Real phi = this->ModifyLimiter( phir, phil );

            qTry[ iEqu ] += phi * ( dqdx * dx + dqdy * dy + dqdz * dz );
        }

        if ( ( * this->ckfun )( qTry ) )
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf2 )[ iEqu ][ ug.fId ] = qTry[ iEqu ];
            }
        }
    }
}

void LimField::CalcFaceValueWeighted()
{
    RealField qTry( this->nEqu );

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        Real dxl = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        Real dyl = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        Real dzl = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        Real dxr = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        Real dyr = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        Real dzr = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        Real delt1  = DIST( dxl, dyl, dzl );
        Real delt2  = DIST( dxr, dyr, dzr );
        Real delta  = 1.0 / ( delt1 + delt2 + SMALL );

        Real cl = delt2 * delta;
        Real cr = delt1 * delta;

        Real dx = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.lc ];
        Real dy = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.lc ];
        Real dz = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.lc ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            qTry[ iEqu ] = ( * this->qf1 )[ iEqu ][ ug.fId ];
        }

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real dqdx = ( * this->dqdx )[ iEqu ][ ug.lc ];
            Real dqdy = ( * this->dqdy )[ iEqu ][ ug.lc ];
            Real dqdz = ( * this->dqdz )[ iEqu ][ ug.lc ];

            Real phil  = ( * this->limiter )[ iEqu ][ ug.lc ];
            Real phir  = ( * this->limiter )[ iEqu ][ ug.rc ];
            Real phi = this->ModifyLimiter( phil, phir );

            qTry[ iEqu ] += phi * ( dqdx * dx + dqdy * dy + dqdz * dz );
        }

        if ( ( * this->ckfun )( qTry ) )
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf1 )[ iEqu ][ ug.fId ] = qTry[ iEqu ];
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                Real var1 = ( * this->qf1 )[ iEqu ][ ug.fId ];
                Real var2 = ( * this->qf2 )[ iEqu ][ ug.fId ];
                Real value = cl * var1 + cr * var2;

                ( * this->qf1 )[ iEqu ][ ug.fId ] = value;
            }
        }

        dx = ( * ug.xfc )[ ug.fId ] - ( * ug.xcc )[ ug.rc ];
        dy = ( * ug.yfc )[ ug.fId ] - ( * ug.ycc )[ ug.rc ];
        dz = ( * ug.zfc )[ ug.fId ] - ( * ug.zcc )[ ug.rc ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            qTry[ iEqu ] = ( * this->qf2 )[ iEqu ][ ug.fId ];
        }

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real dqdx = ( * this->dqdx )[ iEqu ][ ug.rc ];
            Real dqdy = ( * this->dqdy )[ iEqu ][ ug.rc ];
            Real dqdz = ( * this->dqdz )[ iEqu ][ ug.rc ];

            Real phil = ( * this->limiter )[ iEqu ][ ug.lc ];
            Real phir = ( * this->limiter )[ iEqu ][ ug.rc ];
            Real phi = this->ModifyLimiter( phir, phil );

            qTry[ iEqu ] += phi * ( dqdx * dx + dqdy * dy + dqdz * dz );
        }

        if ( ( * this->ckfun )( qTry ) )
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf2 )[ iEqu ][ ug.fId ] = qTry[ iEqu ];
            }
        }
        else
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                Real var1 = ( * this->qf1 )[ iEqu ][ ug.fId ];
                Real var2 = ( * this->qf2 )[ iEqu ][ ug.fId ];
                Real value = cl * var1 + cr * var2;

                ( * this->qf2 )[ iEqu ][ ug.fId ] = value;
            }
        }
    }
}

Limiter::Limiter()
{
    lim = new Lim();
}

Limiter::~Limiter()
{
    delete lim;
}

void Limiter::CalcLimiter()
{
    ug.Init();
    limf->Init();
    Alloc();
    for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
    {
        lim->limiter = & ( * limf->limiter )[ iEqu ];
        lim->q       = & ( * limf->q       )[ iEqu ];
        lim->dqdx    = & ( * limf->dqdx    )[ iEqu ];
        lim->dqdy    = & ( * limf->dqdy    )[ iEqu ];
        lim->dqdz    = & ( * limf->dqdz    )[ iEqu ];
        this->SetInitValue();
        this->CalcLimiterScalar();
    }
    DeAlloc();
}

void Limiter::Alloc()
{
    lim->minvf = new RealField( ug.nTCell );
    lim->maxvf = new RealField( ug.nTCell );
}

void Limiter::DeAlloc()
{
    delete lim->minvf;
    delete lim->maxvf;
}

void Limiter::SetInitValue()
{
    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        ( * lim->limiter )[ cId ] = 1.0;
    }
}

void Limiter::CalcLimiterScalar()
{
    if ( limflag == ILMT_ZERO )
    {
        this->CalcZeroLimiter();
    }
    else if ( limflag == ILMT_NO )
    {
        this->CalcNoLimiter();
    }
    else if ( limflag == ILMT_BARTH )
    {
        this->CalcBarthLimiter();
    }
    else if ( limflag == ILMT_VENCAT )
    {
        this->CalcVencatLimiter();
    }
}

void Limiter::CalcZeroLimiter()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        lim->lim1 = 0.0;
        lim->lim2 = 0.0;

        ( * lim->limiter )[ ug.lc ] = lim->lim1;
        ( * lim->limiter )[ ug.rc ] = lim->lim2;
    }
}

void Limiter::CalcNoLimiter()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        lim->lim1 = 1.0;
        lim->lim2 = 1.0;

        ( * lim->limiter )[ ug.lc ] = lim->lim1;
        ( * lim->limiter )[ ug.rc ] = lim->lim2;
    }
}

void Limiter::CalcBarthLimiter()
{
    this->CalcMinMaxDiff();

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareData();

        this->CalcLocalBarthLimiter();

        ( * lim->limiter )[ ug.lc ] = lim->lim1;
        ( * lim->limiter )[ ug.rc ] = lim->lim2;
    }
}

void Limiter::CalcVencatLimiter()
{
    this->CalcMinMaxDiff();

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareData();

        this->CalcLocalVencatLimiter();

        ( * lim->limiter )[ ug.lc ] = lim->lim1;
        ( * lim->limiter )[ ug.rc ] = lim->lim2;
    }
}

void Limiter::PrepareData()
{
    gcom.xcc1 = ( * ug.xcc )[ ug.lc ];
    gcom.ycc1 = ( * ug.ycc )[ ug.lc ];
    gcom.zcc1 = ( * ug.zcc )[ ug.lc ];

    gcom.xcc2 = ( * ug.xcc )[ ug.rc ];
    gcom.ycc2 = ( * ug.ycc )[ ug.rc ];
    gcom.zcc2 = ( * ug.zcc )[ ug.rc ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    gcom.xfn   = ( * ug.xfn )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn )[ ug.fId ];

    gcom.xfc   = ( * ug.xfc )[ ug.fId ];
    gcom.yfc   = ( * ug.yfc )[ ug.fId ];
    gcom.zfc   = ( * ug.zfc )[ ug.fId ];

    lim->minv1 = ( * lim->minvf )[ ug.lc ];
    lim->minv2 = ( * lim->minvf )[ ug.rc ];

    lim->maxv1 = ( * lim->maxvf )[ ug.lc ];
    lim->maxv2 = ( * lim->maxvf )[ ug.rc ];

    lim->dqdx1 = ( * lim->dqdx )[ ug.lc ];
    lim->dqdy1 = ( * lim->dqdy )[ ug.lc ];
    lim->dqdz1 = ( * lim->dqdz )[ ug.lc ];

    lim->dqdx2 = ( * lim->dqdx )[ ug.rc ];
    lim->dqdy2 = ( * lim->dqdy )[ ug.rc ];
    lim->dqdz2 = ( * lim->dqdz )[ ug.rc ];

    lim->lim1 =  ( * lim->limiter )[ ug.lc ];
    lim->lim2 =  ( * lim->limiter )[ ug.rc ];
}

void Limiter::CalcLocalBarthLimiter()
{
    Real dx1 = gcom.xfc - gcom.xcc1;
    Real dy1 = gcom.yfc - gcom.ycc1;
    Real dz1 = gcom.zfc - gcom.zcc1;

    Real dqFace1  = lim->dqdx1 * dx1 + lim->dqdy1 * dy1 + lim->dqdz1 * dz1;

    Real ds1  = DIST( dx1, dy1, dz1 );
    Real dot1 = ( gcom.xfn * dx1 + gcom.yfn * dy1 + gcom.zfn * dz1 ) / ( ds1 + SMALL );

    Real limv1 = BarthFunction( dqFace1, lim->minv1, lim->maxv1, dot1 );
    lim->lim1 = MIN( lim->lim1, limv1 );

    Real dx2 = gcom.xfc - gcom.xcc2;
    Real dy2 = gcom.yfc - gcom.ycc2;
    Real dz2 = gcom.zfc - gcom.zcc2;

    Real dqFace2  = lim->dqdx2 * dx2 + lim->dqdy2 * dy2 + lim->dqdz2 * dz2;

    Real ds2  = DIST( dx2, dy2, dz2 );
    Real dot2 = - ( gcom.xfn * dx2 + gcom.yfn * dy2 + gcom.zfn * dz2 ) / ( ds2 + SMALL );

    Real limv2 = BarthFunction( dqFace2, lim->minv2, lim->maxv2, dot2 );
    lim->lim2 = MIN( lim->lim2, limv2 );
}

void Limiter::CalcLocalVencatLimiter()
{
    Real eps = ctrl.vencat_coef * MAX( lim->qmax - lim->qmin, 1.0 );
    Real eps1 = SQR( eps ) + SMALL;
    Real eps2 = eps1;

    Real dx1 = gcom.xfc - gcom.xcc1;
    Real dy1 = gcom.yfc - gcom.ycc1;
    Real dz1 = gcom.zfc - gcom.zcc1;

    Real dqFace1  = lim->dqdx1 * dx1 + lim->dqdy1 * dy1 + lim->dqdz1 * dz1;

    Real limv1 = VencatFunction( dqFace1, lim->minv1, lim->maxv1, eps1 );
    lim->lim1 = MIN( lim->lim1, limv1 );

    Real dx2 = gcom.xfc - gcom.xcc2;
    Real dy2 = gcom.yfc - gcom.ycc2;
    Real dz2 = gcom.zfc - gcom.zcc2;

    Real dqFace2  = lim->dqdx2 * dx2 + lim->dqdy2 * dy2 + lim->dqdz2 * dz2;
    Real limv2 = VencatFunction( dqFace2, lim->minv2, lim->maxv2, eps2 );
    lim->lim2 = MIN( lim->lim2, limv2 );
}

void Limiter::CalcMinMaxDiff()
{
    // Find the maximum and minimum in the neighbor of each cell
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * lim->minvf )[ cId ] = ( * lim->q )[ cId ];
        ( * lim->maxvf )[ cId ] = ( * lim->q )[ cId ];
    }

    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        int bcType = ug.bcRecord->bcType[ fId ];
        if ( ! BC::IsInterfaceBc( bcType ) ) continue;

        ( * lim->minvf )[ ug.lc ] = MIN( ( * lim->minvf )[ ug.lc ], ( * lim->q )[ ug.rc ] );
        ( * lim->maxvf )[ ug.lc ] = MAX( ( * lim->maxvf )[ ug.lc ], ( * lim->q )[ ug.rc ] );
    }

    for ( int fId = ug.nBFace; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        ( * lim->minvf )[ ug.lc ] = MIN( ( * lim->minvf )[ ug.lc ], ( * lim->q )[ ug.rc ] );
        ( * lim->maxvf )[ ug.lc ] = MAX( ( * lim->maxvf )[ ug.lc ], ( * lim->q )[ ug.rc ] );

        ( * lim->minvf )[ ug.rc ] = MIN( ( * lim->minvf )[ ug.rc ], ( * lim->q )[ ug.lc ] );
        ( * lim->maxvf )[ ug.rc ] = MAX( ( * lim->maxvf )[ ug.rc ], ( * lim->q )[ ug.lc ] );
    }

    // Get the maximum and the minimum difference
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * lim->minvf )[ cId ] -= ( * lim->q )[ cId ];
        ( * lim->maxvf )[ cId ] -= ( * lim->q )[ cId ];
    }

    lim->qmin =   LARGE;
    lim->qmax = - LARGE;
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        lim->qmin = MIN( lim->qmin, ( * lim->minvf )[ cId ] );
        lim->qmax = MAX( lim->qmax, ( * lim->maxvf )[ cId ] );
    }
}

bool NoCheck( RealField & q )
{
    return true;
}

EndNameSpace