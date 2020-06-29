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

#include "CurveLine.h"
#include "LineMesh.h"
#include "CurveMachine.h"
#include "SegmentCtrl.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

StrCurveLoop::StrCurveLoop()
{
    ;
}

StrCurveLoop::~StrCurveLoop()
{
    ;
}

void StrCurveLoop::AddCurve( int curveId )
{
    this->curveIdList.push_back( curveId );
}

CurveLine * StrCurveLoop::GetCurve( int curveId )
{
    int id = this->curveIdList[ curveId ];
    return curve_Machine.GetCurve( id );
}

void StrCurveLoop::SetDimension()
{
    CurveLine * s1 = this->GetCurve( 0 );
    CurveLine * s2 = this->GetCurve( 1 );
    CurveLine * s3 = this->GetCurve( 2 );
    CurveLine * s4 = this->GetCurve( 3 );

    s1->Alloc( ni );
    s2->Alloc( ni );

    s3->Alloc( nj );
    s4->Alloc( nj );
}

CurveLine::CurveLine()
{
    lineType = 0;
    this->totalCurveLength = 1;
    this->ds1 = -1;
    this->ds2 = -1;
}

CurveLine::~CurveLine()
{
    ;
}

void CurveLine::Alloc( int nNode )
{
    this->nNode = nNode;
    this->x.resize( nNode );
    this->y.resize( nNode );
    this->z.resize( nNode );
}

void CurveLine::MakeLine( PointType & p1, PointType & p2 )
{
    int nPoint = this->x.size();
    for ( int i = 0; i < nPoint; ++ i )
    {
        this->x[ i ] = p1.x + i * ( p2.x - p1.x ) / ( nPoint - 1 );
        this->y[ i ] = p1.y + i * ( p2.y - p1.y ) / ( nPoint - 1 );
        this->z[ i ] = p1.z + i * ( p2.z - p1.z ) / ( nPoint - 1 );
    }
}

void CurveLine::MakeCircle()
{
    PointType & p1 = this->start_p;
    PointType & p2 = this->end_p;
    PointType & p0 = this->center_p;
    this->MakeCircle( p1, p2, p0 );
}

void CurveLine::MakeCircle( PointType & p1, PointType & p2, PointType & p0 )
{
    int nPoint = this->x.size();
    Real dx1 = p1.x - p0.x;
    Real dy1 = p1.y - p0.y;
    Real dz1 = p1.z - p0.z;
    Real r1 = DIST( dx1, dy1, dz1 );

    Real dx2 = p2.x - p0.x;
    Real dy2 = p2.y - p0.y;
    Real dz2 = p2.z - p0.z;

    Real r2 = DIST( dx2, dy2, dz2 );

    Real rad = half * ( r1 + r2 );

    Real div = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
    Real div_unit = div / ( r1 * r2 );

    Real angle_span = acos( div_unit );
    Real dcit = angle_span / ( nPoint - 1 );

    Real cosc1 = dx1 / r1;
    Real sinc1 = dy1 / r1;

    for ( int i = 0; i < nPoint; ++ i )
    {
        //Real dc = i * dcit;
        Real dc = - i * dcit;
        Real sindc = sin( dc );
        Real cosdc = cos( dc );
        Real coscc = cosc1 * cosdc - sinc1 * sindc;
        Real sincc = sinc1 * cosdc + cosc1 * sindc;
        this->x[ i ] = p0.x + rad * coscc;
        this->y[ i ] = p0.y + rad * sincc;
        this->z[ i ] = p0.z;
    }
}

void CurveLine::CalcTotalLength()
{
    if ( lineType == 0 )
    {
        CalcLineLength();
    }
}

void CurveLine::CalcLineLength()
{
    this->totalCurveLength = ONEFLOW::CalcDist( start_p, end_p );
}

void CurveLine::CalcNormal( RealField & nbx, RealField & nby, RealField & nbz )
{
    int nPoint = this->x.size();
    nbx.resize( nPoint );
    nby.resize( nPoint );
    nbz.resize( nPoint );

    Real nnx = 0.0;
    Real nny = 0.0;
    Real nnz = 1.0;

    RealField nref( 3 );
    RealField dtau( 3 );
    RealField norm( 3 );

    nref[ 0 ] = 0.0;
    nref[ 1 ] = 0.0;
    nref[ 2 ] = 1.0;

    for ( int i = 0; i < nPoint - 1; ++ i )
    {
        //tangent;
        Real dx = x[ i + 1 ] - x[ i ];
        Real dy = y[ i + 1 ] - y[ i ];
        Real dz = z[ i + 1 ] - z[ i ];
        Real ods = 1.0 / ( DIST( dx, dy, dz ) );

        dx *= ods;
        dy *= ods;
        dz *= ods;

        dtau[ 0 ] = dx;
        dtau[ 1 ] = dy;
        dtau[ 2 ] = dz;

        CrossProduct( nref, dtau, norm );

        nbx[ i ] = norm[ 0 ];
        nby[ i ] = norm[ 1 ];
        nbz[ i ] = norm[ 2 ];
    }

    nbx[ nPoint - 1 ] = nbx[ nPoint - 2 ];
    nby[ nPoint - 1 ] = nby[ nPoint - 2 ];
    nbz[ nPoint - 1 ] = nbz[ nPoint - 2 ];

    RealField tmpx, tmpy, tmpz;
    tmpx = nbx;
    tmpy = nby;
    tmpz = nbz;

    for ( int i = 0; i < nPoint - 1; ++ i )
    {
        nbx[ i ] = half * ( tmpx[ i ] + tmpx[ i + 1 ] );
        nby[ i ] = half * ( tmpy[ i ] + tmpy[ i + 1 ] );
        nbz[ i ] = half * ( tmpz[ i ] + tmpz[ i + 1 ] );
    }

    int kkk = 1;
}

void CurveLine::GenerateCircleLine()
{
    int nPoint = this->x.size();
    Real dx1 = start_p.x - center_p.x;
    Real dy1 = start_p.y - center_p.y;
    Real dz1 = start_p.z - center_p.z;
    Real r1 = DIST( dx1, dy1, dz1 );

    Real dx2 = end_p.x - center_p.x;
    Real dy2 = end_p.y - center_p.y;
    Real dz2 = end_p.z - center_p.z;

    Real r2 = DIST( dx2, dy2, dz2 );

    Real rad = half * ( r1 + r2 );

    Real div = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
    Real div_unit = div / ( r1 * r2 );

    Real angle_span = acos( div_unit );
    Real dcit = angle_span / ( nPoint - 1 );

    Real cosc1 = dx1 / r1;
    Real sinc1 = dy1 / r1;

    Real ave_ratio = 1.0 / ( nPoint - 1 );

    Real ratio = 0.6 * ave_ratio;

    Real coef = 1.0;
    Real cc = 1.0 / ( nPoint - 1 );
    ONEFLOW::GetExponentialCoeff( ratio, cc, coef );

    for ( int i = 0; i < nPoint; ++ i )
    {
        Real r = static_cast< Real >( i ) / ( nPoint - 1 );
        Real ar = coef * r;
        Real factor = ( exp( ar ) - 1 ) / ( exp( coef ) - 1 );

        Real dc = - angle_span * factor;
        Real sindc = sin( dc );
        Real cosdc = cos( dc );
        Real coscc = cosc1 * cosdc - sinc1 * sindc;
        Real sincc = sinc1 * cosdc + cosc1 * sindc;

        this->x[ i ] = center_p.x + rad * coscc;
        this->y[ i ] = center_p.y + rad * sincc;
        this->z[ i ] = center_p.z;
    }
}

Real CurveLine::CalcParabolicLength(Real x, Real p)
{
    Real term = 2 * x / p;
    Real s = p / 2 * (sqrt(term * (1 + term)) + log(sqrt(term) + sqrt(1 + term)));
    return s;
}

void CurveLine::FindXByLength(Real & x, Real p, Real maxX, Real length)
{
    Real x1 = 0;
    Real x2 = maxX;

    while (true)
    {
        x = half * (x1 + x2);

        Real s = this->CalcParabolicLength(x, p);

        if (s < length)
        {
            x1 = x;
        }
        else
        {
            x2 = x;
        }

        if (ABS(s - length) < 1.0e-8 || ABS(x2 - x1) < 1.0e-8)
        {
            break;
        }
    };
}

void CurveLine::GenerateParabolicLine()
{
    int nPoint = this->x.size();
    Real L = ABS( start_p.x - end_p.x );
    Real H = ABS( start_p.y - end_p.y );
    Real p = H * H / (2.0 * L);

    Real curve_lenth = this->CalcParabolicLength( L, p );

    Real ave_ratio = 1.0 / (nPoint - 1);

    Real ratio = 0.5 * ave_ratio;

    Real coef = 1.0;
    Real cc = 1.0 / ( nPoint - 1 );
    ONEFLOW::GetExponentialCoeff(ratio, cc, coef);

    for (int i = 0; i < nPoint; ++i)
    {
        Real r = static_cast<Real>(i) / (nPoint - 1);
        Real ar = coef * r;
        Real factor = (exp(ar) - 1) / (exp(coef) - 1);

        Real local_len = curve_lenth * factor;

        Real xlocal;
        this->FindXByLength( xlocal, p, L, local_len );
        Real ylocal = sqrt(2 * p * xlocal);

        this->x[i] = xlocal - L;
        this->y[i] = ylocal;
        this->z[i] = center_p.z;
    }
}

void CurveLine::GenerateCurveLine()
{
    CalcLineLength();

    int nPoint = this->x.size();
    this->ds1 = this->totalCurveLength / ( nPoint - 1 );
    Real ratio = ds1 / this->totalCurveLength;
    Real coef = 0.0001;
    Real cc = 1.0 / ( nPoint - 1 );
    ONEFLOW::GetExponentialCoeff( ratio, cc, coef );

    Real dx = end_p.x - start_p.x;
    Real dy = end_p.y - start_p.y;
    Real dz = end_p.z - start_p.z;
    for ( int i = 0; i < nPoint; ++ i )
    {
        Real r = static_cast< Real >( i ) / ( nPoint - 1 );
        Real ar = coef * r;
        Real factor = ( exp( ar ) - 1 ) / ( exp( coef ) - 1 );

        Real xx = start_p.x + factor * dx;
        Real yy = start_p.y + factor * dy;
        Real zz = start_p.z + factor * dz;

        this->x[ i ] = xx;
        this->y[ i ] = yy;
        this->z[ i ] = zz;
    }
}

Real CalcDist( PointType & p1, PointType & p2 )
{
    Real dx = p1.x - p2.x;
    Real dy = p1.y - p2.y;
    Real dz = p1.z - p2.z;
    Real ds = DIST( dx, dy, dz );
    return ds;
}

void CrossProduct( RealField & a, RealField & b, RealField & c )
{
    c[ 0 ] = a[ 1 ] * b[ 2 ] - a[ 2 ] * b[ 1 ];
    c[ 1 ] = a[ 2 ] * b[ 0 ] - a[ 0 ] * b[ 2 ];
    c[ 2 ] = a[ 0 ] * b[ 1 ] - a[ 1 ] * b[ 0 ];
}


EndNameSpace