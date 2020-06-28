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

#include "SegmentCtrl.h"
#include "LineMachine.h"
#include "LineMesh.h"
#include "FileIO.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

SegmentCtrl::SegmentCtrl()
{
    this->ratio1 = -1;
    this->ratio2 = -1;
    segmentCopy = 0;
}

SegmentCtrl::~SegmentCtrl()
{
    delete segmentCopy;
}

void SegmentCtrl::Read( FileIO * ioFile )
{
    string distributionString = ioFile->ReadNextWord();

    if ( distributionString == "r" )
    {
        this->distribution = 0;
        this->ratio1 = ioFile->ReadNextDigit< Real >();
        this->ratio2 = ioFile->ReadNextDigit< Real >();
    }
    else if ( distributionString == "d" )
    {
        this->distribution = 1;
        this->ds1 = ioFile->ReadNextDigit< Real >();
        this->ds2 = ioFile->ReadNextDigit< Real >();
    }
    else if ( distributionString.substr( 0, 1 ) == "c" )
    {
        this->distribution = 2;
        this->segmentCopy = new SegmentCopy();
        this->segmentCopy->Read( ioFile );
    }
    else if ( distributionString.substr( 0, 1 ) == "e" )
    {
        this->distribution = 3;
        this->cA1 = 0.5;
        this->cA2 = 1.0e-4;
        this->cA3 = 0.5;
    }
}


void SegmentCtrl::CalcFactor()
{
    if ( this->distribution == 2 )
    {
        this->CopyFactor();
    }
    else
    {
        this->CalcEffectiveFactor();
    }
}

Real SegmentCtrl::CalcFactor( Real compCoor )
{
    Real factor = -1;
    if ( this->distribution == 3 )
    {
        factor = this->CalcDFactor( compCoor );
    }
    else
    {
        factor = this->CalcDFactor( compCoor );
        //factor = this->CalcSFactor( compCoor );

    }
    return factor;
}

Real SegmentCtrl::CalcDFactor( Real compCoor )
{
    Real factor;
    if ( compCoor <= cA3 )
    {
        Real term1 = ( cA2 / cA3 ) * compCoor;
        Real factor1 = cA1 * ( exp( term1 ) - 1.0 ) / ( exp( cA2 ) - 1.0 );
        factor = factor1;
    }
    else
    {
        Real term2 = cA4 * ( compCoor - cA3 ) / ( 1.0 - cA3 );
        Real factor2 = cA1 + ( 1.0 - cA1 ) * ( exp( term2 ) - 1.0 ) / ( exp( cA4 ) - 1.0 );
        factor = factor2;
    }
    return factor;
}

Real SegmentCtrl::CalcSFactor( Real compCoor )
{
    Real cFactor = ( exp( expCoeff * compCoor ) - 1.0 ) / ( exp( expCoeff ) - 1.0 );
    Real factor = c1 * cFactor + c2 * ( 1 - cFactor );
    return factor;
}

void SegmentCtrl::InitCoef()
{
    ;
}

void SegmentCtrl::InitDExp()
{
    Real term1 = cA1 * cA2 / cA3;
    Real term2 = ( 1.0 - cA3 ) / ( 1.0 - cA1 );
    Real term3 = exp( cA2 ) / ( exp( cA2 ) - 1.0 );
    Real term = term1 * term2 * term3;

    Real xsta, xend;

    xsta = - 20.0;
    xend = 20.0;
    Real fmid;
    Real xmid;
    while ( true )
    {
        xmid = half * ( xsta + xend );
        if ( xmid == 0 )
        {
            fmid = 1.0;
        }
        else
        {
            fmid = xmid / ( exp( xmid ) - 1.0 );
        }

        if ( ABS( fmid - term ) < 1.0e-8 ) break;

        if ( fmid < term )
        {
            xend = xmid;
        }
        else
        {
            xsta = xmid;
        }
    };

    cA4 = xmid;
}

Real SegmentCtrl::CalDExp()
{
    Real term1 = cA1 * cA2 / cA3;
    Real term2 = ( 1.0 - cA3 ) / ( 1.0 - cA1 );
    Real term3 = exp( cA2 ) / ( exp( cA2 ) - 1.0 );
    Real term = term1 * term2 * term3;
    Real fmid = cA4 / ( exp( cA4 ) - 1.0 );
    Real diff = fmid - term;
    return diff;
}

void SegmentCtrl::Init()
{
    //this->CalcEffectiveRatio();
    this->CalcEffectiveRatioTest();

    //if ( this->distribution == 3 )
    //{
    //    this->InitDExp();
    //}
    //else
    //{
    //    Real ratio_new = this->ratio / ( this->nPoint - 1 );
    //    Real cc = 1.0 / ( nPoint - 1 );
    //    GetExponentialCoeff( ratio_new, cc, expCoeff );
    //}
}

void SegmentCtrl::SetPara( Real diff, Real & v_min, Real & v_max )
{
    //if ( diff > 0.0 )
    //{
    //    v_min = cA1;
    //}
    //else
    //{
    //    v_max = cA1;
    //}
    //this->cA1 = half * ( v_min + v_max );
    //cout << " cA1 = " << cA1 << " diff = " << diff << "\n";

    if ( diff < 0.0 )
    {
        v_min = cA3;
    }
    else
    {
        v_max = cA3;
    }
    this->cA3 = half * ( v_min + v_max );
    cout << " cA3 = " << cA3 << " diff = " << diff << "\n";
}

void SegmentCtrl::CalcEffectiveRatioTest()
{
    if ( ds1 < 0.0 && ds2 < 0.0 )
    {
        this->c1 = 1;
        this->cA1 = 1.0;
        this->cA2 = 1.0;
        this->cA3 = 1.0;
        this->cA4 = 1.0;

        //Real rds1 = 1.0 / ( this->nPoint - 1 );
        //Real fds1 = rds1 / cA1;
        //Real rho1 = 1.0 / ( this->nPoint - 1 );
        //Real cc1 = rho1 / cA3;
        //GetExponentialCoeff( fds1, cc1, cA2 );
    }
    else if ( ds1 > 0.0 && ds2 < 0.0 )
    {
        this->c1 = 1;
        this->cA1 = 1.0;
        this->cA3 = 1.0;

        Real rds1 = ds1 / lenth;
        Real fds1 = rds1 / cA1;
        Real rho1 = 1.0 / ( this->nPoint - 1 );
        Real cc1 = rho1 / cA3;
        GetExponentialCoeff( fds1, cc1, cA2 );
    }
    else if ( ds1 < 0.0 && ds2 > 0.0 )
    {
        this->c1 = 0;
        this->cA1 = 1.0;
        this->cA3 = 1.0;
        Real rds1 = ds2 / lenth;
        Real fds1 = rds1 / cA1;
        Real rho1 = 1.0 / ( this->nPoint - 1 );
        Real cc1 = rho1 / cA3;
        GetExponentialCoeff( fds1, cc1, cA2 );
    }
    else if ( ds1 > 0.0 && ds2 > 0.0 )
    {
        this->c1 = 1;
        this->cA1 = 0.5;
        this->cA3 = 0.5;
        Real v_min = 0.001;
        Real v_max = 0.999;
        this->cA1 = half * ( v_min + v_max );
        int icount = 0;
        while ( true )
        {
            Real rds1 = ds1 / lenth;
            Real rds2 = ds2 / lenth;
            Real fds1 = rds1 / cA1;
            Real fds2 = 1 - rds2 / ( 1 - cA1 );
            Real rho1 = 1.0 / ( this->nPoint - 1 );
            Real rho2 = ( this->nPoint - 2.0 ) / ( this->nPoint - 1 );
            Real cc1 = rho1 / cA3;
            Real cc2 = ( rho2 - cA3 ) / ( 1 - cA3 );

            GetExponentialCoeff( fds1, cc1, cA2 );
            GetExponentialCoeff( fds2, cc2, cA4 );

            Real diff = CalDExp();
            ++ icount;
            if ( ABS( diff ) < 1.0e-4 || icount > 100 ) break;

            SetPara( diff, v_min, v_max );
        
        }

    }

    this->c2 = 1 - this->c1;
}

void SegmentCtrl::CalcEffectiveRatio()
{
    if ( this->distribution == 1 )
    {
        Real ds = this->lenth / ( this->nPoint - 1 );
        this->ratio1 = this->ds1 / ds;
        this->ratio2 = this->ds2 / ds;
    }

    this->ratio = 1.0;

    this->c1 = 1;

    if ( ratio1 > 0.0 && ratio2 < 0.0 )
    {
        this->ratio = this->ratio1;
    }
    else if ( ratio1 <= 0.0 && ratio2 > 0.0 )
    {
        this->ratio = this->ratio2;
        this->c1 = 0;
    }

    this->c2 = 1 - this->c1;
}

void SegmentCtrl::CalcEffectiveFactor()
{
    this->Init();

    this->factorList.resize( nPoint );
    this->pidxList.resize( nPoint );

    for ( int iPoint = 1; iPoint < nPoint - 1; ++ iPoint )
    {
        Real coor = static_cast<Real>( iPoint ) / ( nPoint - 1 );
        Real factor = this->CalcFactor( coor );

        int idx = c1 * iPoint + c2 * ( nPoint - 1 - iPoint );
        this->factorList[ iPoint ] = factor;
        this->pidxList[ iPoint ] = idx;
    }
}

void SegmentCtrl::CopyFactor()
{
    this->factorList.resize( nPoint );
    this->pidxList.resize( nPoint );
    int nTPoint = segmentCopy->GetNPoint();
    RealField dsArray( nTPoint );
    Real totalS = 0.0;
    for ( int iPoint = 1; iPoint < nTPoint; ++ iPoint )
    {
        PointType & p0 = segmentCopy->GetPoint( iPoint - 1 );
        PointType & p1 = segmentCopy->GetPoint( iPoint );
        Real x0 = p0.x;
        Real y0 = p0.y;
        Real z0 = p0.z;

        Real x1 = p1.x;
        Real y1 = p1.y;
        Real z1 = p1.z;

        Real dx = x1 - x0;
        Real dy = y1 - y0;
        Real dz = z1 - z0;
        Real ds = DIST( dx, dy, dz );
        totalS += ds;
        dsArray[ iPoint ] = totalS;
    }

    for ( int iPoint = 1; iPoint < nPoint - 1; ++ iPoint )
    {
        Real factor = dsArray[ iPoint ] / totalS;
        int effectivePointIndex = iPoint;
        this->factorList[ iPoint ] = factor;
        this->pidxList[ iPoint ] = iPoint;
    }
}

SegmentCopy::SegmentCopy()
{
    ;
}

SegmentCopy::~SegmentCopy()
{
    ;
}

void SegmentCopy::Read( FileIO * ioFile )
{
    while ( true )
    {
        string word = ioFile->ReadNextWord();
        if ( word == "" ) break;
        int line_id = StringToDigit< int >( word );
        lineList.push_back( line_id );
    }
}

int SegmentCopy::GetNPoint()
{
    int nSegment = this->GetNSegment();
    this->lowerIdList.resize( nSegment );
    this->upperIdList.resize( nSegment );

    int lowerIndex = 0;
    int upperIndex = 0;

    int nTPoint = 0;
    for ( int iSegment = 0; iSegment < nSegment; ++ iSegment )
    {
        int id = lineList[ iSegment ];
        CurveMesh * CurveMesh = line_Machine.GetCurveMesh( id );
        int nPoint = CurveMesh->segmentCtrl->nPoint;
        lowerIndex = upperIndex;
        upperIndex = lowerIndex + nPoint - 1;

        this->lowerIdList[ iSegment ] = lowerIndex;
        this->upperIdList[ iSegment ] = upperIndex;

        nTPoint = upperIndex + 1;
    }
    return nTPoint;
}

PointType & SegmentCopy::GetPoint( int pointIndex, int signFlag )
{
    int segmentId = -1;
    int localId = -1;
    this->FindSegmentAndLocalId( pointIndex, signFlag, segmentId, localId );
    int id = lineList[ segmentId ];
    int ida = ABS( id );
    CurveMesh * curveMesh = line_Machine.GetCurveMesh( ida );
    int localSignFlag = 1;
    if ( id < 0 ) localSignFlag = -1;

    return curveMesh->GetPoint( localId, localSignFlag );
}

void SegmentCopy::FindSegmentAndLocalId( int id, int signFlag, int & segmentId, int & localId )
{
    int nPoint = this->GetNPoint();
    int coef = ( 1 + signFlag ) / 2;
    int fid = coef * id + ( 1 - coef ) * ( nPoint - 1 - id );

    int nSegment = this->GetNSegment();
    for ( int iSegment = 0; iSegment < nSegment; ++ iSegment )
    {
        if ( this->lowerIdList[ iSegment ] <= fid &&
            fid <= this->upperIdList[ iSegment ] )
        {
            segmentId = iSegment;
            localId = fid - this->lowerIdList[ iSegment ];
            return;
        }
    }
}

void GetExponentialCoeff( Real targetRatio, Real cc, Real & coef )
{
    Real coefMin = -100.0;
    Real coefMax = 100.0;

    Real ratio;

    int iCount = 0;
    while ( true )
    {
        coef = half * ( coefMin + coefMax );
        if ( ABS( coef ) < 1.0e-8 )
        {
            ratio = cc;
        }
        else
        {
            ratio = ( exp( coef * cc ) - 1.0 ) / ( exp( coef ) - 1.0 );
        }
        
        if ( ABS( ( ratio - targetRatio ) ) / targetRatio < 1.0e-4 || iCount > 100 )
        {
            break;
        }
        if ( ratio > targetRatio )
        {
            coefMin = coef;
        }
        else
        {
            coefMax = coef;
        }
        ++ iCount;
    }
}


EndNameSpace