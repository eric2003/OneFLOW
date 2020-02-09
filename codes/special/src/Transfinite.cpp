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

#include "Transfinite.h"
#include "CurveLine.h"
#include "CurveMesh.h"
#include "DataBaseIO.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void TransfiniteInterpolation( RealField2D & coor, int ni, int nj )
{
    RealField2D u, v, uv;

    AllocateVector( u , ni, nj );
    AllocateVector( v , ni, nj );
    AllocateVector( uv, ni, nj );

    int il = 0;
    int ir = ni - 1;

    int jl = 0;
    int jr = nj - 1;

    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            Real u0 = coor[ il ][ j ];
            Real u1 = coor[ ir ][ j ];

            Real v0 = coor[ i ][ jl ];
            Real v1 = coor[ i ][ jr ];

            Real uv00 = coor[ il ][ jl ];
            Real uv01 = coor[ il ][ jr ];
            Real uv10 = coor[ ir ][ jl ];
            Real uv11 = coor[ ir ][ jr ];

            Real coefu = static_cast< Real >( i ) / ( ni - 1 );
            Real coefv = static_cast< Real >( j ) / ( nj - 1 );

            Real coefuv00 = ( 1 - coefu ) * ( 1 - coefv );
            Real coefuv01 = ( 1 - coefu ) * ( coefv     );
            Real coefuv10 = ( coefu     ) * ( 1 - coefv );
            Real coefuv11 = ( coefu     ) * ( coefv     );

            u [ i ][ j ] = ( 1 - coefu ) * u0 + coefu * u1;
            v [ i ][ j ] = ( 1 - coefv ) * v0 + coefv * v1;
            uv[ i ][ j ] = coefuv00 * uv00 +
                           coefuv01 * uv01 +
                           coefuv10 * uv10 +
                           coefuv11 * uv11;
        }
    }

    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            Real uu   = u [ i ][ j ];
            Real vv   = v [ i ][ j ];
            Real uuvv = uv[ i ][ j ];
            coor[ i ][ j ] = uu + vv - uuvv;
        }
    }
}

void TransfiniteInterpolation( RealField3D & coor, int ni, int nj, int nk )
{
    RealField3D u, v, w;
    RealField3D uv, uw, vw;
    RealField3D uvw;

    AllocateVector( u, ni, nj, nk );
    AllocateVector( v, ni, nj, nk );
    AllocateVector( w, ni, nj, nk );
    AllocateVector( uv, ni, nj, nk );
    AllocateVector( uw, ni, nj, nk );
    AllocateVector( vw, ni, nj, nk );
    AllocateVector( uvw, ni, nj, nk );

    int il = 0;
    int ir = ni - 1;

    int jl = 0;
    int jr = nj - 1;

    int kl = 0;
    int kr = nk - 1;
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            for ( int i = 0; i < ni; ++ i )
            {
                Real u0 = coor[ il ][ j ][ k ];
                Real u1 = coor[ ir ][ j ][ k ];

                Real v0 = coor[ i ][ jl ][ k ];
                Real v1 = coor[ i ][ jr ][ k ];

                Real w0 = coor[ i ][ j ][ kl ];
                Real w1 = coor[ i ][ j ][ kr ];

                Real x0j0 = coor[ il ][ j ][ kl ];
                Real x0j1 = coor[ il ][ j ][ kr ];
                Real x1j0 = coor[ ir ][ j ][ kl ];
                Real x1j1 = coor[ ir ][ j ][ kr ];

                Real x00k = coor[ il ][ jl ][ k ];
                Real x01k = coor[ il ][ jr ][ k ];
                Real x10k = coor[ ir ][ jl ][ k ];
                Real x11k = coor[ ir ][ jr ][ k ];

                Real xi00 = coor[ i ][ jl ][ kl ];
                Real xi01 = coor[ i ][ jl ][ kr ];
                Real xi10 = coor[ i ][ jr ][ kl ];
                Real xi11 = coor[ i ][ jr ][ kr ];

                Real x000 = coor[ il ][ jl ][ kl ];
                Real x100 = coor[ ir ][ jl ][ kl ];
                Real x110 = coor[ ir ][ jr ][ kl ];
                Real x010 = coor[ il ][ jr ][ kl ];

                Real x001 = coor[ il ][ jl ][ kr ];
                Real x101 = coor[ ir ][ jl ][ kr ];
                Real x111 = coor[ ir ][ jr ][ kr ];
                Real x011 = coor[ il ][ jr ][ kr ];

                Real ci = static_cast<Real>( i ) / ( ni - 1 );
                Real cj = static_cast<Real>( j ) / ( nj - 1 );
                Real ck = static_cast<Real>( k ) / ( nk - 1 );

                Real di = 1 - ci;
                Real dj = 1 - cj;
                Real dk = 1 - ck;

                u[ i ][ j ][ k ] = di * u0 + ci * u1;
                v[ i ][ j ][ k ] = dj * v0 + cj * v1;
                w[ i ][ j ][ k ] = dk * v0 + ck * w1;

                uw[ i ][ j ][ k ] = di * dk * x0j0 + di * ck * x0j1 + ci * dk * x1j0 + ci * ck * x1j1;
                uv[ i ][ j ][ k ] = di * dj * x00k + di * cj * x01k + ci * dj * x10k + ci * ck * x11k;
                vw[ i ][ j ][ k ] = dj * dk * xi00 + dj * ck * xi01 + cj * dk * xi10 + cj * ck * xi11;
                uvw[ i ][ j ][ k ] = di * dj * dk * x000 + ci * dj * dk * x100 + ci * cj * dk * x110 + di * cj * dk * x010
                                   + di * dj * ck * x001 + ci * dj * ck * x101 + ci * cj * ck * x111 + di * cj * ck * x011;

            }
        }
    }
    for ( int k = 1; k < nk - 1; ++ k )
    {
        for ( int j = 1; j < nj - 1; ++ j )
        {
            for ( int i = 1; i < ni - 1; ++ i )
            {
                Real u0 = u[ i ][ j ][ k ];
                Real v0 = v[ i ][ j ][ k ];
                Real w0 = w[ i ][ j ][ k ];
                Real uv0 = uv[ i ][ j ][ k ];
                Real vw0 = vw[ i ][ j ][ k ];
                Real uw0 = uw[ i ][ j ][ k ];
                Real uvw0 = uvw[ i ][ j ][ k ];

                coor[ i ][ j ][ k ] = u0 + v0 + w0 - uv0 - vw0 - uw0 + uvw0;
            }
        }
    }
}

void AlgebraInterpolation( RealField2D & x, RealField2D & y, RealField2D & z, int ni, int nj, 
    RealField & nbx, RealField & nby, RealField & nbz, Real beta )
{
    //for ( int i = 1; i < ni - 1; ++ i )
    for ( int i = 0; i < ni; ++ i )
    {
        Real x1 = x[ i ][ 0 ];
        Real y1 = y[ i ][ 0 ];
        Real z1 = z[ i ][ 0 ];

        Real xN = x[ i ][ nj - 1 ];
        Real yN = y[ i ][ nj - 1 ];
        Real zN = z[ i ][ nj - 1 ];

        Real dx = xN - x1;
        Real dy = yN - y1;
        Real dz = zN - z1;
        Real ds = DIST( dx, dy, dz );

        Real nx = nbx[ i ];
        Real ny = nby[ i ];
        Real nz = nbz[ i ];

        Real rbx = dx / ds;
        Real rby = dy / ds;
        Real rbz = dz / ds;

        for ( int j = 0; j < nj; ++ j )
        {
            Real b = static_cast< Real >( j ) /( nj - 1 );
            Real term1 = ( beta + 1 ) / ( beta - 1 );
            Real term2 = pow( term1, 1 - b );
            Real ak = 1 + beta * ( 1 - term2 ) / ( 1 + term2 );
            Real fk = 1 - ak;

            Real cx = ak * ( nx * fk + rbx *( 1 - fk ) );
            Real cy = ak * ( ny * fk + rby *( 1 - fk ) );
            Real cz = ak * ( nz * fk + rbz *( 1 - fk ) );
            x[ i ][ j ] = x1 + ds * cx;
            y[ i ][ j ] = y1 + ds * cy;
            z[ i ][ j ] = z1 + ds * cz;
        }
    }
}


EndNameSpace