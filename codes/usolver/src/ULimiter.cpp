/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "UGradPolicy.h"
#include "ULimiterPolicy.h"
#include "ULimiterRange.h"
#include "HipProductionPolicy.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "Zone.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "Ctrl.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "Iteration.h"

#include <algorithm>

BeginNameSpace( ONEFLOW )

void LimField::BindFieldsAndAllocFaces(
    UnsGrid * grid,
    const char * qName,
    const char * dqdxName,
    const char * dqdyName,
    const char * dqdzName,
    const char * limiterName )
{
    q       = GetFieldPointer< MRField >( grid, qName );
    dqdx    = GetFieldPointer< MRField >( grid, dqdxName );
    dqdy    = GetFieldPointer< MRField >( grid, dqdyName );
    dqdz    = GetFieldPointer< MRField >( grid, dqdzName );
    limiter = GetFieldPointer< MRField >( grid, limiterName );

    this->nEqu = q->GetNEqu();
    // Init() runs once per RK stage.  qf1/qf2 are fully overwritten by
    // GetQlQr()/CalcFaceValue() before every read, so reuse them whenever the
    // extent still matches instead of reallocating ~9 MB per stage.
    const HXSize_t nFaces = grid->nFaces;
    const bool needsAllocation =
        qf1 == nullptr || qf2 == nullptr
        || qf1->GetNEqu() != this->nEqu
        || qf2->GetNEqu() != this->nEqu
        || this->nEqu <= 0
        || ( * qf1 )[ 0 ].size() != nFaces
        || ( * qf2 )[ 0 ].size() != nFaces;
    if ( needsAllocation )
    {
        qf1 = std::make_unique< MRField >( this->nEqu, nFaces );
        qf2 = std::make_unique< MRField >( this->nEqu, nFaces );
    }
}


Real LimField::ModifyLimiter( Real phil, Real phir )
{
    return BlendFaceLimiter( phil, phir, ug.ireconface != 0 );
}

void LimField::GetQlQr()
{
    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            ( * this->qf1 )[ iEqu ][ ug.fId ] = ( * this->q )[ iEqu ][ ug.lc ];
            ( * this->qf2 )[ iEqu ][ ug.fId ] = ( * this->q )[ iEqu ][ ug.rc ];
        }
    }
}

void LimField::BcQlQrFix()
{
    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        int bcType = ug.bcRecord->bcType[ fId ];
        if ( bcType == BC::INTERFACE ) continue;
        if ( bcType == BC::PERIODIC  ) continue;

        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real tmp = ArithmeticFaceAverage(
                ( * this->q )[ iEqu ][ ug.lc ], ( * this->q )[ iEqu ][ ug.rc ] );

            ( * this->qf1 )[ iEqu ][ ug.fId ] = tmp;
            ( * this->qf2 )[ iEqu ][ ug.fId ] = tmp;
        }
    }
}


void LimField::CalcFaceValue()
{
    RealField qTry( this->nEqu );

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

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

            qTry[ iEqu ] = ApplyLimitedFaceReconstruction(
                qTry[ iEqu ], phi, dqdx, dqdy, dqdz, dx, dy, dz );
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

            qTry[ iEqu ] = ApplyLimitedFaceReconstruction(
                qTry[ iEqu ], phi, dqdx, dqdy, dqdz, dx, dy, dz );
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

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
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

        Real cl, cr;
        InverseDistanceFaceWeights( dxl, dyl, dzl, dxr, dyr, dzr, cl, cr );

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

            qTry[ iEqu ] = ApplyLimitedFaceReconstruction(
                qTry[ iEqu ], phi, dqdx, dqdy, dqdz, dx, dy, dz );
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
                Real value = WeightedFaceAverage( var1, var2, cl, cr );

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

            qTry[ iEqu ] = ApplyLimitedFaceReconstruction(
                qTry[ iEqu ], phi, dqdx, dqdy, dqdz, dx, dy, dz );
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
                Real value = WeightedFaceAverage( var1, var2, cl, cr );

                ( * this->qf2 )[ iEqu ][ ug.fId ] = value;
            }
        }
    }
}

void Limiter::CalcLimiter()
{
    ug.Init();
    limf->Init();

    if ( limflag == ILMT_NO )
    {
        // No-limiter mode must leave phi == 1 for every cell and equation, and
        // CalcLimiterScalar() performs no face traversal in that mode anyway.
        // Fill the field directly instead of allocating and releasing the
        // min/max scratch pair once per equation; the values are identical
        // (exactly one, no arithmetic), so the CPU oracle is unchanged.
        // Repeating the fill every stage was the only per-step host loop left in
        // the residual chain, but it also repaired the ghosts that the interface
        // exchange overwrites.  Skip the repeat only when that exchange no
        // longer carries the field.
        if ( noLimiterFieldFilled ) return;
        for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
        {
            RealField & limiterField = ( * limf->limiter )[ iEqu ];
            std::fill( limiterField.begin(), limiterField.end(), one );
        }
        noLimiterFieldFilled = HipUnusedInterfaceFieldsSkipped();
        return;
    }

    // Keep the extrema workspace across calls. assign() clears every entry,
    // including ghost-cell slots, while reusing the vector capacity when the
    // grid extent is unchanged. The zero initialization matches the previous
    // freshly-constructed vectors and is required by the current boundary path.
    if ( ! lim->minvf ) lim->minvf = std::make_unique<RealField>();
    if ( ! lim->maxvf ) lim->maxvf = std::make_unique<RealField>();
    lim->minvf->assign( ug.nTCell, 0.0 );
    lim->maxvf->assign( ug.nTCell, 0.0 );

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
    // The owned extrema workspace remains available for the next limiter call.
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
    // SetInitValue() has already filled the complete cell field with one.
    // The no-limiter mode therefore needs no face traversal at all.
    if ( limflag == ILMT_NO ) return;

    if ( limflag == ILMT_ZERO )
    {
        this->CalcZeroLimiter();
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
    for ( int fId = 0; fId < ug.nFaces; ++ fId )
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

void Limiter::CalcBarthLimiter()
{
    this->BuildNeighborRangeOffsets();

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
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
    this->BuildNeighborRangeOffsets();

    // The global offset range is fixed for this equation during the face loop.
    const Real regularizationEpsilon = VencatRegularizationEpsilon(
        ctrl.vencat_coef, lim->minOffsetRange, lim->maxOffsetRange );

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareData();

        this->CalcLocalVencatLimiter( regularizationEpsilon );

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
    const Real dx1 = gcom.xfc - gcom.xcc1;
    const Real dy1 = gcom.yfc - gcom.ycc1;
    const Real dz1 = gcom.zfc - gcom.zcc1;

    const Real dqFace1 = ReconstructFaceIncrement(
        lim->dqdx1, lim->dqdy1, lim->dqdz1, dx1, dy1, dz1 );
    const Real dot1 = FaceSideDirectionDot(
        gcom.xfn, gcom.yfn, gcom.zfn, dx1, dy1, dz1, 1.0 );
    lim->lim1 = TightenLimiter(
        lim->lim1, BarthFunction( dqFace1, lim->minv1, lim->maxv1, dot1 ) );

    const Real dx2 = gcom.xfc - gcom.xcc2;
    const Real dy2 = gcom.yfc - gcom.ycc2;
    const Real dz2 = gcom.zfc - gcom.zcc2;

    const Real dqFace2 = ReconstructFaceIncrement(
        lim->dqdx2, lim->dqdy2, lim->dqdz2, dx2, dy2, dz2 );
    const Real dot2 = FaceSideDirectionDot(
        gcom.xfn, gcom.yfn, gcom.zfn, dx2, dy2, dz2, -1.0 );
    lim->lim2 = TightenLimiter(
        lim->lim2, BarthFunction( dqFace2, lim->minv2, lim->maxv2, dot2 ) );
}

void Limiter::CalcLocalVencatLimiter( Real regularizationEpsilon )
{
    const Real dx1 = gcom.xfc - gcom.xcc1;
    const Real dy1 = gcom.yfc - gcom.ycc1;
    const Real dz1 = gcom.zfc - gcom.zcc1;

    const Real dqFace1 = ReconstructFaceIncrement(
        lim->dqdx1, lim->dqdy1, lim->dqdz1, dx1, dy1, dz1 );
    lim->lim1 = TightenLimiter(
        lim->lim1,
        VencatFunction( dqFace1, lim->minv1, lim->maxv1, regularizationEpsilon ) );

    const Real dx2 = gcom.xfc - gcom.xcc2;
    const Real dy2 = gcom.yfc - gcom.ycc2;
    const Real dz2 = gcom.zfc - gcom.zcc2;

    const Real dqFace2 = ReconstructFaceIncrement(
        lim->dqdx2, lim->dqdy2, lim->dqdz2, dx2, dy2, dz2 );
    lim->lim2 = TightenLimiter(
        lim->lim2,
        VencatFunction( dqFace2, lim->minv2, lim->maxv2, regularizationEpsilon ) );
}

void Limiter::BuildNeighborRangeOffsets()
{
    BuildCellNeighborOffsetBounds(
        * lim->q,
        * ug.lcf,
        * ug.rcf,
        ug.bcRecord->bcType,
        ug.nCells,
        ug.nBFaces,
        ug.nFaces,
        []( int bcType ) { return BC::IsInterfaceBc( bcType ); },
        * lim->minvf,
        * lim->maxvf,
        lim->minOffsetRange,
        lim->maxOffsetRange );
}

MRField * Limiter::GetLeftField() const
{
    return limf->qf1.get();
}

MRField * Limiter::GetRightField() const
{
    return limf->qf2.get();
}

void Limiter::CalcFaceValue()
{
    limf->CalcFaceValue();
}

void Limiter::GetQlQr()
{
    limf->GetQlQr();
}

void Limiter::BcQlQrFix()
{
    limf->BcQlQrFix();
}


bool NoCheck( RealField & q )
{
    return true;
}

EndNameSpace
