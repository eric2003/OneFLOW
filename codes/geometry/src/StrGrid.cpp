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

#include "StrGrid.h"
#include "Dimension.h"
#include "BcRecord.h"
#include "InterFace.h"
#include "HXMath.h"
#include "HXMathExt.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
#include "DataBaseIO.h"
#include "DataBook.h"
#include "Tolerence.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_GRID( StrGrid )

StrGrid::StrGrid()
{
    faceTopo = new FaceTopo();
    bcRegionGroup = new BcRegionGroup();
    strx = 0;
    stry = 0;
    strz = 0;
}

StrGrid::~StrGrid()
{
    delete faceTopo;
    delete bcRegionGroup;
    delete strx;
    delete stry;
    delete strz;
}

void StrGrid::Decode( DataBook * databook )
{
    this->ReadGrid( databook );
}

void StrGrid::Encode( DataBook * databook )
{
    this->WriteGrid( databook );
}

int StrGrid::CalcNumberOfNode()
{
    if ( ONEFLOW::IsTwoD() )
    {
        return ONEFLOW::CalcNumberOfNode( ni, nj );
    }
    return ONEFLOW::CalcNumberOfNode( ni, nj, nk );
}

int StrGrid::CalcNumberOfCell()
{
    if ( ONEFLOW::IsTwoD() )
    {
        return ONEFLOW::CalcNumberOfCell( ni, nj );
    }
    return ONEFLOW::CalcNumberOfCell( ni, nj, nk );
}

int StrGrid::CalcNumberOfFace()
{
    if ( ONEFLOW::IsTwoD() )
    {
        return ONEFLOW::CalcNumberOfFace( ni, nj );
    }
    return ONEFLOW::CalcNumberOfFace( ni, nj, nk );
}

void StrGrid::SetBasicDimension()
{
    this->nNode = this->CalcNumberOfNode();
    this->nCell = this->CalcNumberOfCell();
    this->nFace = this->CalcNumberOfFace();

    cout << " " << Dim::dimension << "D Grid\n";

    cout << " number of nodes           : " << this->nNode << endl;
    cout << " number of element surface : " << this->nFace << endl;
    cout << " number of element         : " << this->nCell << endl;
}

void StrGrid::SetLayout()
{
   if ( this->nodeMesh->xN.size() != 0 &&
        this->nodeMesh->yN.size() != 0 &&
        this->nodeMesh->zN.size() != 0 )
    {
        Range I, J, K;
        I.SetRange( 1, ni );
        J.SetRange( 1, nj );
        K.SetRange( 1, nk );

        this->strx = new Field3D( & this->nodeMesh->xN[ 0 ], I, J, K );
        this->stry = new Field3D( & this->nodeMesh->yN[ 0 ], I, J, K );
        this->strz = new Field3D( & this->nodeMesh->zN[ 0 ], I, J, K );
    }
}

void StrGrid::ReadGrid( DataBook * databook )
{
}

void StrGrid::ReadBoundaryTopology( DataBook * databook )
{
}

void StrGrid::WriteGrid( DataBook * databook )
{
}

void StrGrid::WriteGridFaceTopology( DataBook * databook )
{
}

void StrGrid::WriteBoundaryTopology( DataBook * databook )
{
}

void StrGrid::CalcMinMaxDis3D( Real & dismin, Real & dismax )
{
    dismin =   LARGE;
    dismax = - LARGE;

    Field3D & xs = * this->strx;
    Field3D & ys = * this->stry;
    Field3D & zs = * this->strz;

    int ist, ied, jst, jed, kst, ked;

    ist = 1;
    ied = ni - 1;

    jst = 1;
    jed = nj - 1;

    kst = 1;
    ked = nk - 1;
    const int nEdge = 12;
    Real dx[ nEdge ], dy[ nEdge ], dz[ nEdge ];

    Real ptTol = Tolerence::GetTol();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                dx[ 0 ]  = xs( i + 1, j  , k  ) - xs( i  , j  , k  );
                dy[ 0 ]  = ys( i + 1, j  , k  ) - ys( i  , j  , k  );
                dz[ 0 ]  = zs( i + 1, j  , k  ) - zs( i  , j  , k  );

                dx[ 1 ]  = xs( i + 1, j + 1, k  ) - xs( i + 1, j  , k  );
                dy[ 1 ]  = ys( i + 1, j + 1, k  ) - ys( i + 1, j  , k  );
                dz[ 1 ]  = zs( i + 1, j + 1, k  ) - zs( i + 1, j  , k  );

                dx[ 2 ]  = xs( i  , j + 1, k  ) - xs( i + 1, j + 1, k  );
                dy[ 2 ]  = ys( i  , j + 1, k  ) - ys( i + 1, j + 1, k  );
                dz[ 2 ]  = zs( i  , j + 1, k  ) - zs( i + 1, j + 1, k  );

                dx[ 3 ]  = xs( i  , j  , k  ) - xs( i  , j + 1, k  );
                dy[ 3 ]  = ys( i  , j  , k  ) - ys( i  , j + 1, k  );
                dz[ 3 ]  = zs( i  , j  , k  ) - zs( i  , j + 1, k  );

                dx[ 4 ]  = xs( i + 1, j  , k + 1 ) - xs( i  , j  , k + 1 );
                dy[ 4 ]  = ys( i + 1, j  , k + 1 ) - ys( i  , j  , k + 1 );
                dz[ 4 ]  = zs( i + 1, j  , k + 1 ) - zs( i  , j  , k + 1 );

                dx[ 5 ]  = xs( i + 1, j + 1, k + 1 ) - xs( i + 1, j  , k + 1 );
                dy[ 5 ]  = ys( i + 1, j + 1, k + 1 ) - ys( i + 1, j  , k + 1 );
                dz[ 5 ]  = zs( i + 1, j + 1, k + 1 ) - zs( i + 1, j  , k + 1 );

                dx[ 6 ]  = xs( i    , j + 1, k + 1 ) - xs( i + 1, j + 1, k + 1 );
                dy[ 6 ]  = ys( i    , j + 1, k + 1 ) - ys( i + 1, j + 1, k + 1 );
                dz[ 6 ]  = zs( i    , j + 1, k + 1 ) - zs( i + 1, j + 1, k + 1 );

                dx[ 7 ]  = xs( i    , j    , k + 1 ) - xs( i    , j + 1, k + 1 );
                dy[ 7 ]  = ys( i    , j    , k + 1 ) - ys( i    , j + 1, k + 1 );
                dz[ 7 ]  = zs( i    , j    , k + 1 ) - zs( i    , j + 1, k + 1 );

                dx[ 8 ]  = xs( i    , j    , k + 1 ) - xs( i    , j    , k    );
                dy[ 8 ]  = ys( i    , j    , k + 1 ) - ys( i    , j    , k    );
                dz[ 8 ]  = zs( i    , j    , k + 1 ) - zs( i    , j    , k    );

                dx[ 9 ]  = xs( i + 1, j    , k + 1 ) - xs( i + 1, j    , k    );
                dy[ 9 ]  = ys( i + 1, j    , k + 1 ) - ys( i + 1, j    , k    );
                dz[ 9 ]  = zs( i + 1, j    , k + 1 ) - zs( i + 1, j    , k    );

                dx[ 10 ] = xs( i + 1, j + 1, k + 1 ) - xs( i + 1, j + 1, k    );
                dy[ 10 ] = ys( i + 1, j + 1, k + 1 ) - ys( i + 1, j + 1, k    );
                dz[ 10 ] = zs( i + 1, j + 1, k + 1 ) - zs( i + 1, j + 1, k    );

                dx[ 11 ] = xs( i    , j + 1, k + 1 ) - xs( i    , j + 1, k    );
                dy[ 11 ] = ys( i    , j + 1, k + 1 ) - ys( i    , j + 1, k    );
                dz[ 11 ] = zs( i    , j + 1, k + 1 ) - zs( i    , j + 1, k    );

                for ( int m = 0; m < nEdge; ++ m )
                {
                    Real ds = ONEFLOW::DIST( dx[ m ], dy[ m ], dz[ m ] );
                    if ( ds <= ptTol ) continue;
                    dismin = ONEFLOW::MIN( dismin, ds );
                    dismax = ONEFLOW::MAX( dismax, ds );
                }
            }
        }
    }
}

void StrGrid::CalcMinMaxDis2D( Real & dismin, Real & dismax )
{
    dismin =   LARGE;
    dismax = - LARGE;

    Field3D & xs = * this->strx;
    Field3D & ys = * this->stry;
    Field3D & zs = * this->strz;

    int ist, ied, jst, jed, kst, ked;

    ist = 1;
    ied = ni - 1;

    jst = 1;
    jed = nj - 1;

    kst = 1;
    ked = 1;
    const int nEdge = 4;
    Real dx[ nEdge ], dy[ nEdge ], dz[ nEdge ];

    Real ptTol = Tolerence::GetTol();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                dx[ 0 ]  = xs( i + 1, j    , k ) - xs( i    , j    , k );
                dy[ 0 ]  = ys( i + 1, j    , k ) - ys( i    , j    , k );
                dz[ 0 ]  = zs( i + 1, j    , k ) - zs( i    , j    , k );

                dx[ 1 ]  = xs( i + 1, j + 1, k ) - xs( i + 1, j    , k );
                dy[ 1 ]  = ys( i + 1, j + 1, k ) - ys( i + 1, j    , k );
                dz[ 1 ]  = zs( i + 1, j + 1, k ) - zs( i + 1, j    , k );

                dx[ 2 ]  = xs( i    , j + 1, k ) - xs( i + 1, j + 1, k );
                dy[ 2 ]  = ys( i    , j + 1, k ) - ys( i + 1, j + 1, k );
                dz[ 2 ]  = zs( i    , j + 1, k ) - zs( i + 1, j + 1, k );

                dx[ 3 ]  = xs( i    , j    , k ) - xs( i    , j + 1, k );
                dy[ 3 ]  = ys( i    , j    , k ) - ys( i    , j + 1, k );
                dz[ 3 ]  = zs( i    , j    , k ) - zs( i    , j + 1, k );

                for ( int m = 0; m < nEdge; ++ m )
                {
                    Real ds = ONEFLOW::DIST( dx[ m ], dy[ m ], dz[ m ] );
                    if ( ds <= ptTol ) continue;
                    dismin = ONEFLOW::MIN( dismin, ds );
                    dismax = ONEFLOW::MAX( dismax, ds );
                }
            }
        }
    }
}

void StrGrid::CalcMinMaxDis1D( Real & dismin, Real & dismax )
{
    dismin =   LARGE;
    dismax = - LARGE;

    Field3D & xs = * this->strx;
    Field3D & ys = * this->stry;
    Field3D & zs = * this->strz;

    int ist, ied, jst, jed, kst, ked;

    ist = 1;
    ied = ni - 1;

    jst = 1;
    jed = 1;

    kst = 1;
    ked = 1;
    const int nEdge = 1;
    Real dx[ nEdge ], dy[ nEdge ], dz[ nEdge ];

    Real ptTol = Tolerence::GetTol();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                dx[ 0 ]  = xs( i + 1, j    , k ) - xs( i    , j    , k );
                dy[ 0 ]  = ys( i + 1, j    , k ) - ys( i    , j    , k );
                dz[ 0 ]  = zs( i + 1, j    , k ) - zs( i    , j    , k );

                for ( int m = 0; m < nEdge; ++ m )
                {
                    Real ds = ONEFLOW::DIST( dx[ m ], dy[ m ], dz[ m ] );
                    if ( ds <= ptTol ) continue;
                    dismin = ONEFLOW::MIN( dismin, ds );
                    dismax = ONEFLOW::MAX( dismax, ds );
                }
            }
        }
    }
}

void StrGrid::GetMinMaxDistance( Real & dismin, Real & dismax )
{
    if ( Dim::dimension == THREE_D )
    {
        this->CalcMinMaxDis3D( dismin, dismax );
    }
    else if ( Dim::dimension == TWO_D )
    {
        this->CalcMinMaxDis2D( dismin, dismax );
    }
    else
    {
        this->CalcMinMaxDis1D( dismin, dismax );
    }
}

int CalcNumberOfFace( const int & ni )
{
    return ni;
}

int CalcNumberOfFace( const int & ni, const int & nj )
{
    return ( nj - 1 ) * ni + ( ni - 1 ) * nj;
}

int CalcNumberOfFace( const int & ni, const int & nj, const int & nk )
{
    return ( nj - 1 ) * ( nk - 1 ) * ni + ( nk - 1 ) * ( ni - 1 ) * nj + ( ni - 1 ) * ( nj - 1 ) * nk;
}

int CalcNumberOfNode( const int & ni )
{
    return ni;
}

int CalcNumberOfNode( const int & ni, const int & nj )
{
    return ni * nj;
}

int CalcNumberOfNode( const int & ni, const int & nj, const int & nk )
{
    return ni * nj * nk;
}

int CalcNumberOfCell( const int & ni )
{
    return ONEFLOW::COUNT( 1, ni - 1 );
}

int CalcNumberOfCell( const int & ni, const int & nj )
{
    return ONEFLOW::CalcNumberOfCell( ni ) * ONEFLOW::CalcNumberOfCell( nj );
}

int CalcNumberOfCell( const int & ni, const int & nj, const int & nk )
{
    return ONEFLOW::CalcNumberOfCell( ni ) *
           ONEFLOW::CalcNumberOfCell( nj ) *
           ONEFLOW::CalcNumberOfCell( nk ) ;
}

StrGrid * StrGridCast( Grid * gridIn )
{
    return static_cast< StrGrid * >( gridIn );
}

void GetIJKRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K )
{
    I.SetRange( 1 + startShift, ni + endShift );
    J.SetRange( 1 + startShift, nj + endShift );
    K.SetRange( 1 + startShift, nk + endShift );

    if ( ni == 1 ) I.SetRange( 1, 1 );
    if ( nj == 1 ) J.SetRange( 1, 1 );
    if ( nk == 1 ) K.SetRange( 1, 1 );
}
Range IJKRange::I;
Range IJKRange::J;
Range IJKRange::K;

int IJKRange::ist;
int IJKRange::ied;
int IJKRange::jst;
int IJKRange::jed;
int IJKRange::kst;
int IJKRange::ked;

IJKRange::IJKRange()
{
    ;
}

IJKRange::~IJKRange()
{
    ;
}

void IJKRange::Calc( int ni, int nj, int nk, int ss, int es )
{
    IJKRange::I.SetRange( 1 + ss, ni + es );
    IJKRange::J.SetRange( 1 + ss, nj + es );
    IJKRange::K.SetRange( 1 + ss, nk + es );

    if ( ni == 1 ) IJKRange::I.SetRange( 1, 1 );
    if ( nj == 1 ) IJKRange::J.SetRange( 1, 1 );
    if ( nk == 1 ) IJKRange::K.SetRange( 1, 1 );
}

void IJKRange::ToScalar()
{
    IJKRange::ist = IJKRange::I.First();
    IJKRange::ied = IJKRange::I.Last();
    IJKRange::jst = IJKRange::J.First();
    IJKRange::jed = IJKRange::J.Last();
    IJKRange::kst = IJKRange::K.First();
    IJKRange::ked = IJKRange::K.Last();
}

EndNameSpace