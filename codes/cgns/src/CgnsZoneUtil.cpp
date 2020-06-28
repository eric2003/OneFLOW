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

#include "CgnsZone.h"
#include "CgnsZoneUtil.h"
#include "CgnsBase.h"
#include "CgnsCoor.h"
#include "CgnsSection.h"
#include "CgnsZsection.h"
#include "CgnsBcBoco.h"
#include "CgnsZbc.h"
#include "CgnsZbcBoco.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "NodeMesh.h"
#include "StrUtil.h"
#include "Grid.h"
#include "BgGrid.h"
#include "StrGrid.h"
#include "GridState.h"
#include "Dimension.h"
#include "GridElem.h"
#include "ElemFeature.h"
#include "ElementHome.h"
#include "PointFactory.h"
#include "PointSearch.h"
#include "FaceSolver.h"
#include "Stop.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

void EncodeIJK( int & index, int i, int j, int k, int ni, int nj, int nk )
{
    index = ( i - 1 ) + ( j - 1 ) * ni + ( k - 1 ) * ( ni * nj ) ;
}

void DecodeIJK( int index, int & i, int & j, int & k, int ni, int nj, int nk )
{
    k = index / ( ni * nj ) + 1;
    index -= ( k - 1 ) * ni * nj;
    j = index / ni + 1;
    i = index - ( j - 1 ) * ni + 1;
}

void GetRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K )
{
    I.SetRange( 1 + startShift, ni + endShift );
    J.SetRange( 1 + startShift, nj + endShift );
    K.SetRange( 1 + startShift, nk + endShift );

    if ( ni == 1 ) I.SetRange( 1, 1 );
    if ( nj == 1 ) J.SetRange( 1, 1 );
    if ( nk == 1 ) K.SetRange( 1, 1 );
}

void GetIJKRegion( Range & I, Range & J, Range & K, int & ist, int & ied, int & jst, int & jed, int & kst, int & ked )
{
    ist = I.First();
    ied = I.Last();
    jst = J.First();
    jed = J.Last();
    kst = K.First();
    ked = K.Last();
}

void PrepareCgnsZoneSub( Grids & grids, CgnsZone * cgnsZone )
{
    NodeMesh * nodeMesh = cgnsZone->cgnsCoor->GetNodeMesh();

    int nNode, nCell;

    HXVector< Int3D * > unsIdList;

    MergeToSingleZone( grids, unsIdList, nodeMesh, nNode, nCell );

    cgnsZone->cgnsCoor->SetNNode( nNode );
    cgnsZone->cgnsCoor->SetNCell( nCell );

    FillSection( grids, unsIdList, cgnsZone );

    cgnsZone->ConvertToInnerDataStandard();

    ONEFLOW::DeletePointer( unsIdList );
}

void MergeToSingleZone( Grids & grids, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell )
{
    PointSearch * point_search = new PointSearch();
    point_search->Initialize( grids );

    size_t nZone = grids.size();

    unsIdList.resize( nZone );
    nCell = 0;
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;
        nCell += grid->nCell;
        unsIdList[ iZone ] = new Int3D( Range( 1, ni ), Range( 1, nj ), Range( 1, nk ) );
        Int3D & unsId = * unsIdList[ iZone ];
        cout << " block = " << iZone + 1 << "\n";
        CalcUnsId( grid, point_search, & unsId );
    }

    nNode = point_search->GetNPoint();

    cout << " First nNode = " << nNode << "\n";
    nodeMesh->xN.resize( nNode );
    nodeMesh->yN.resize( nNode );
    nodeMesh->zN.resize( nNode );
    for ( int i = 0; i < nNode; ++ i )
    {
        Real xm, ym, zm;
        point_search->GetPoint( i, xm, ym, zm );

        nodeMesh->xN[ i ] = xm;
        nodeMesh->yN[ i ] = ym;
        nodeMesh->zN[ i ] = zm;
    }
    delete point_search;
}

void FillSection( Grids & grids, HXVector< Int3D * > & unsIdList, CgnsZone * cgnsZone )
{
    int nTBcRegion = 0;

    int nTCell = 0;
    int nBFace = 0;

    for ( int iZone = 0; iZone < grids.size(); ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        Int3D & unsId = * unsIdList[ iZone ];

        nTCell += grid->CalcNumberOfCell();

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        size_t nBcRegions = bcRegionGroup->regions->size();

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = ( * bcRegionGroup->regions )[ ir ];
            if ( BC::IsNotNormalBc( bcRegion->bcType ) ) continue;
            
            nBFace += bcRegion->CalcRegionCells();
            nTBcRegion ++;
        }
    }

    cout << " nBFace = " << nBFace << "\n";

    cgnsZone->cgnsCoor->SetNCell( nTCell );

    cgnsZone->cgnsZsection->nSection = 2;
    cgnsZone->cgnsZsection->CreateCgnsSection();

    cgnsZone->cgnsZsection->cgnsSections[ 0 ]->startId = 1;
    cgnsZone->cgnsZsection->cgnsSections[ 0 ]->endId   = nTCell;

    cgnsZone->cgnsZsection->cgnsSections[ 1 ]->startId = nTCell + 1;
    cgnsZone->cgnsZsection->cgnsSections[ 1 ]->endId   = nTCell + 1 + nBFace;

    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        cgnsZone->cgnsZsection->cgnsSections[ 0 ]->eType = HEXA_8;
        cgnsZone->cgnsZsection->cgnsSections[ 1 ]->eType = QUAD_4;
    }
    else
    {
        cgnsZone->cgnsZsection->cgnsSections[ 0 ]->eType = QUAD_4;
        cgnsZone->cgnsZsection->cgnsSections[ 1 ]->eType = BAR_2;
    }

    cgnsZone->cgnsZsection->CreateConnList();

    CgnsZbc * cgnsZbc = cgnsZone->cgnsZbc;
    cgnsZbc->cgnsZbcBoco->ReadZnboco( nTBcRegion );
    cgnsZbc->CreateCgnsZbc( cgnsZbc );

    CgnsSection * secV = cgnsZone->cgnsZsection->GetCgnsSection( 0 );
    CgnsSection * secB = cgnsZone->cgnsZsection->GetCgnsSection( 1 );

    CgIntField& connList  = secV->connList;
    CgIntField& bConnList = secB->connList;

    int pos = 0;

    for ( int iZone = 0; iZone < grids.size(); ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        Int3D & unsId = * unsIdList[ iZone ];
        
        IJKRange::Calc( ni, nj, nk, 0, -1 );
        IJKRange::ToScalar();

        int is = 1;
        int js = 1;
        int ks = 1;

        if ( Dim::dimension == ONEFLOW::TWO_D ) ks = 0;

        int eNodeNumbers = ONEFLOW::GetElementNodeNumbers( secV->eType );

        for ( int k = IJKRange::kst; k <= IJKRange::ked; ++ k )
        {
            for ( int j = IJKRange::jst; j <= IJKRange::jed; ++ j )
            {
                for ( int i = IJKRange::ist; i <= IJKRange::ied; ++ i )
                {
                    connList[ pos + 0 ] = unsId( i   , j   , k    ) + 1;
                    connList[ pos + 1 ] = unsId( i+is, j   , k    ) + 1;
                    connList[ pos + 2 ] = unsId( i+is, j+js, k    ) + 1;
                    connList[ pos + 3 ] = unsId( i   , j+js, k    ) + 1;
                    if ( Dim::dimension == ONEFLOW::THREE_D )
                    {
                        connList[ pos + 4 ] = unsId( i   , j   , k+ks ) + 1;
                        connList[ pos + 5 ] = unsId( i+is, j   , k+ks ) + 1;
                        connList[ pos + 6 ] = unsId( i+is, j+js, k+ks ) + 1;
                        connList[ pos + 7 ] = unsId( i   , j+js, k+ks ) + 1;
                    }
                    pos += eNodeNumbers;
                }
            }
        }
    }

    secV->SetElemPosition();
    secB->SetElemPosition();

    int irc  = 0;
    int eIdPos  = nTCell;
    pos = 0;

    BcTypeMap * bcTypeMap = new BcTypeMap();
    bcTypeMap->Init();

    for ( int iZone = 0; iZone < grids.size(); ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        Int3D & unsId = * unsIdList[ iZone ];

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        size_t nBcRegions = bcRegionGroup->regions->size();

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = ( * bcRegionGroup->regions )[ ir ];
            if ( BC::IsNotNormalBc( bcRegion->bcType ) ) continue;
            int nRegionCell = bcRegion->CalcRegionCells();

            CgnsBcBoco * cgnsBcBoco = cgnsZbc->cgnsZbcBoco->GetCgnsBc( irc );
            
            cgnsBcBoco->SetCgnsBcRegionGridLocation( CellCenter );
            cgnsBcBoco->nElements    = 2;
            cgnsBcBoco->bcType       = static_cast< BCType_t >( bcTypeMap->OneFlow2Cgns( bcRegion->bcType ) );
            cgnsBcBoco->pointSetType = PointRange;

            //cgnsBcBoco->SetCgnsBcRegion( nElements, bcType, );

            cgnsBcBoco->CreateCgnsBcConn();
            cgnsBcBoco->connList[ 0 ] = eIdPos + 1;
            cgnsBcBoco->connList[ 1 ] = eIdPos + nRegionCell;
            string bcName = GetCgnsBcName( cgnsBcBoco->bcType );
            cgnsBcBoco->name = AddString( bcName, ir );

            eIdPos += nRegionCell;

            ONEFLOW::SetUnsBcConn( bcRegion, bConnList, pos, unsId );

            irc ++;
        }
    }

    delete bcTypeMap;
}

void CalcUnsId( StrGrid * grid, PointSearch * pointSearch, Int3D * unsId )
{
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;

    Field3D & xs = * grid->strx;
    Field3D & ys = * grid->stry;
    Field3D & zs = * grid->strz;

    RealField coordinate( 3 );
    for ( int k = 1; k <= nk; ++ k )
    {
        for ( int j = 1; j <= nj; ++ j )
        {
            for ( int i = 1; i <= ni; ++ i )
            {
                Real xm = xs( i, j, k );
                Real ym = ys( i, j, k );
                Real zm = zs( i, j, k );

                int pointIndex = pointSearch->AddPoint( xm, ym, zm );

                //{
                //    int width = 8;
                //    cout << " id = " << pointIndex;
                //    cout << setw( width ) << xm;
                //    cout << setw( width ) << ym;
                //    cout << setw( width ) << zm;
                //    cout << "\n";
                //}
                

                ( * unsId )( i, j, k ) = pointIndex;
            }
        }
    }
}

void SetUnsBcConn( BcRegion * bcRegion, CgIntField& conn, int & pos, Int3D & unsId )
{
    int ist, ied, jst, jed, kst, ked;
    bcRegion->GetNormalizeIJKRegion( ist, ied, jst, jed, kst, ked );

    cout << " ist, ied, jst, jed, kst, ked = " << ist << " " << ied << " " << jst << " " << jed << " " << kst << " " << ked << endl;
    int numpt = 4;
    if ( Dim::dimension == TWO_D ) numpt = 2;

    if ( ist == ied )
    {
        int i = ist;
        if ( Dim::dimension == THREE_D )
        {
            for ( int k = kst; k <= ked - 1; ++ k )
            {
                for ( int j = jst; j <= jed - 1; ++ j )
                {
                    if ( i == 1 )
                    {
                        conn[ pos + 0 ] = unsId( i, j    , k     ) + 1;
                        conn[ pos + 1 ] = unsId( i, j    , k + 1 ) + 1;
                        conn[ pos + 2 ] = unsId( i, j + 1, k + 1 ) + 1;
                        conn[ pos + 3 ] = unsId( i, j + 1, k     ) + 1;
                    }
                    else
                    {
                        conn[ pos + 0 ] = unsId( i, j    , k     ) + 1;
                        conn[ pos + 1 ] = unsId( i, j + 1, k     ) + 1;
                        conn[ pos + 2 ] = unsId( i, j + 1, k + 1 ) + 1;
                        conn[ pos + 3 ] = unsId( i, j    , k + 1 ) + 1;
                    }
                    pos += numpt;
                }
            }
        }
        else
        {
            int k = kst;
            for ( int j = jst; j <= jed - 1; ++ j )
            {
                if ( i == 1 )
                {
                    conn[ pos + 0 ] = unsId( i, j + 1, k  ) + 1;
                    conn[ pos + 1 ] = unsId( i, j    , k  ) + 1;
                }
                else
                {
                    conn[ pos + 0 ] = unsId( i, j    , k  ) + 1;
                    conn[ pos + 1 ] = unsId( i, j + 1, k  ) + 1;
                }
                pos += numpt;
            }
        }
        return;
    }

    if ( jst == jed )
    {
        int j = jst;
        if ( Dim::dimension == THREE_D )
        {
            for ( int k = kst; k <= ked - 1; ++ k )
            {
                for ( int i = ist; i <= ied - 1; ++ i )
                {
                    if ( j == 1 )
                    {
                        conn[ pos + 0 ] = unsId( i    , j, k    ) + 1;
                        conn[ pos + 1 ] = unsId( i + 1, j, k    ) + 1;
                        conn[ pos + 2 ] = unsId( i + 1, j, k + 1 ) + 1;
                        conn[ pos + 3 ] = unsId( i    , j, k + 1 ) + 1;
                    }   
                    else
                    {
                        conn[ pos + 0 ] = unsId( i    , j, k     ) + 1;
                        conn[ pos + 1 ] = unsId( i    , j, k + 1 ) + 1;
                        conn[ pos + 2 ] = unsId( i + 1, j, k + 1 ) + 1;
                        conn[ pos + 3 ] = unsId( i + 1, j, k     ) + 1;
                    }
                    pos += numpt;
                }
            }
        }
        else
        {
            int k = kst;
            for ( int i = ist; i <= ied - 1; ++ i )
            {
                if ( j == 1 )
                {
                    conn[ pos + 0 ] = unsId( i    , j, k  ) + 1;
                    conn[ pos + 1 ] = unsId( i + 1, j, k  ) + 1;
                }   
                else
                {
                    conn[ pos + 0 ] = unsId( i + 1, j, k  ) + 1;
                    conn[ pos + 1 ] = unsId( i    , j, k  ) + 1;
                }
                pos += numpt;
            }
        }
        return;
    }

    if ( kst == ked )
    {
        int k = kst;
        for ( int j = jst; j <= jed - 1; ++ j )
        {
            for ( int i = ist; i <= ied - 1; ++ i )
            {
                if ( k == 1 )
                {
                    conn[ pos + 0 ] = unsId( i    , j    , k ) + 1;
                    conn[ pos + 1 ] = unsId( i    , j + 1, k ) + 1;
                    conn[ pos + 2 ] = unsId( i + 1, j + 1, k ) + 1;
                    conn[ pos + 3 ] = unsId( i + 1, j    , k ) + 1;
                }   
                else
                {
                    conn[ pos + 0 ] = unsId( i    , j    , k ) + 1;
                    conn[ pos + 1 ] = unsId( i + 1, j    , k ) + 1;
                    conn[ pos + 2 ] = unsId( i + 1, j + 1, k ) + 1;
                    conn[ pos + 3 ] = unsId( i    , j + 1, k ) + 1;
                }
                pos += numpt;
            }
        }
        return;
    }

    Stop( " error : ist != ied, jst != jed, kst != ked \n" );
}

void GenerateUnsBcElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    int iSection = 1;
    CgnsSection * cgnsSection = myZone->cgnsZsection->GetCgnsSection( iSection );

    myZone->cgnsZbc->CreateCgnsZbc( cgnsZoneIn->cgnsZbc );

    cout << " ConnectionList Size = " << cgnsSection->connSize << "\n";
    cgnsZoneIn->cgnsZbc->GenerateUnsBcElemConn( cgnsSection->connList );
}

void GenerateUnsBcCondConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    int iSection = 1;
    CgnsSection * cgnsSection = myZone->cgnsZsection->GetCgnsSection( iSection );

    CgInt startId = cgnsSection->startId;

    int nBoco = cgnsZoneIn->cgnsZbc->cgnsZbcBoco->nBoco;
    for ( int iBoco = 0; iBoco < nBoco; ++ iBoco )
    {
        CgnsBcBoco * bcRegion    = myZone    ->cgnsZbc->cgnsZbcBoco->GetCgnsBc( iBoco );
        CgnsBcBoco * strBcRegion = cgnsZoneIn->cgnsZbc->cgnsZbcBoco->GetCgnsBc( iBoco );
        bcRegion->CopyStrBcRegion( strBcRegion, startId );
    }
}

void GenerateUnsVolElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    int ni = static_cast<int> (cgnsZoneIn->GetNI());
    int nj = static_cast<int> (cgnsZoneIn->GetNJ());
    int nk = static_cast<int> (cgnsZoneIn->GetNK());

    cout << " ni = " << ni << " nj = " << nj << " nk = " << nk << "\n";

    int iSection = 0;
    CgnsSection * cgnsSection = myZone->cgnsZsection->GetCgnsSection( iSection );

    Range I, J, K;
    GetRange( ni, nj, nk, 0, -1, I, J, K );

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion( I, J, K, ist, ied, jst, jed, kst, ked );

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;

    int cell_dim = myZone->cgnsBase->celldim;

    if ( cell_dim == TWO_D ) kl1 = 0;
    if ( cell_dim == ONE_D ) jl1 = 0;

    CgIntField & connList = cgnsSection->connList;

    int pos = 0;

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                int index1, index2, index3, index4;
                EncodeIJK( index1,  i      , j      , k,  ni,  nj,  nk );
                EncodeIJK( index2,  i + il1, j      , k,  ni,  nj,  nk );

                connList[ pos ++ ] = myZone->l2g[ index1 ] + 1;
                connList[ pos ++ ] = myZone->l2g[ index2 ] + 1;

                if ( cell_dim == ONE_D ) continue;

                EncodeIJK( index3,  i + il1, j + jl1, k,  ni,  nj,  nk );
                EncodeIJK( index4,  i      , j + jl1, k,  ni,  nj,  nk );

                connList[ pos ++ ] = myZone->l2g[ index3 ] + 1;
                connList[ pos ++ ] = myZone->l2g[ index4 ] + 1;

                if ( cell_dim == TWO_D ) continue;

                int index5, index6, index7, index8;
                EncodeIJK( index5,  i      , j      , k + kl1,  ni,  nj,  nk );
                EncodeIJK( index6,  i + il1, j      , k + kl1,  ni,  nj,  nk );
                EncodeIJK( index7,  i + il1, j + jl1, k + kl1,  ni,  nj,  nk );
                EncodeIJK( index8,  i      , j + jl1, k + kl1,  ni,  nj,  nk );

                connList[ pos ++ ] = myZone->l2g[ index5 ] + 1;
                connList[ pos ++ ] = myZone->l2g[ index6 ] + 1;
                connList[ pos ++ ] = myZone->l2g[ index7 ] + 1;
                connList[ pos ++ ] = myZone->l2g[ index8 ] + 1;
            }
        }
    }
}

void AllocateUnsElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    myZone->cgnsZsection->nSection = 2;
    myZone->cgnsZsection->CreateCgnsSection();

    int s1, e1, s2, e2, etype1, etype2;
    //cgnsZoneIn->GetStrZonePara( s1, e1, s2, e2, etype1, etype2 );
    ONEFLOW::GetStrZonePara( cgnsZoneIn, s1, e1, s2, e2, etype1, etype2 );

    CgnsSection * cgnsSection1 = myZone->cgnsZsection->GetCgnsSection( 0 );
    CgnsSection * cgnsSection2 = myZone->cgnsZsection->GetCgnsSection( 1 );
    cgnsSection1->SetSectionInfo( "Section1", etype1, s1, e1 );
    cgnsSection2->SetSectionInfo( "Section2", etype2, s2, e2 );

    myZone->cgnsZsection->CreateConnList();
}

void ReadElementConnectivities( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    ONEFLOW::AllocateUnsElemConn( myZone, cgnsZoneIn );
    ONEFLOW::GenerateUnsVolElemConn( myZone, cgnsZoneIn );
    ONEFLOW::GenerateUnsBcElemConn( myZone, cgnsZoneIn );
    myZone->SetElemPosition();
    ONEFLOW::GenerateUnsBcCondConn( myZone, cgnsZoneIn );
}

void GetStrZonePara( CgnsZone * myZone, int & s1, int & e1, int & s2, int & e2, int & etype1, int & etype2  )
{
    int nActualBcFace = myZone->cgnsZbc->GetNumberOfActualBcElements();

    s1 = 1;
    e1 = myZone->cgnsCoor->GetNCell();

    s2 = e1 + 1;
    e2 = e1 + nActualBcFace;

    int celldim = myZone->cgnsBase->celldim;

    if ( celldim == ONE_D )
    {
        etype1  = CGNS_ENUMV( BAR_2 );
        etype2  = CGNS_ENUMV( NODE );
    }
    else if ( celldim == TWO_D )
    {
        etype1  = CGNS_ENUMV( QUAD_4 );
        etype2  = CGNS_ENUMV( BAR_2  );
    }
    else if ( celldim == THREE_D )
    {
        etype1  = CGNS_ENUMV( HEXA_8 );
        etype2  = CGNS_ENUMV( QUAD_4 );
    }
}

void ReadCgnsZoneType( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    myZone->cgnsZoneType = CGNS_ENUMV( Unstructured );

    cout << "   The Zone Type is " << GetCgnsZoneTypeName( myZone->cgnsZoneType ) << " Zone" << "\n";
}

void ReadCgnsZoneNameAndGeneralizedDimension( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    myZone->zoneName = cgnsZoneIn->zoneName;
}

void SetDimension( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    CgnsCoor * cgnsCoorIn = cgnsZoneIn->cgnsCoor;
    myZone->cgnsCoor->SetDimension( cgnsCoorIn );
}

void ReadCgnsGridCoordinates( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    NodeMesh * nodeMesh1 = myZone->cgnsCoor->GetNodeMesh();
    NodeMesh * nodeMesh2 = cgnsZoneIn->cgnsCoor->GetNodeMesh();

    * nodeMesh1 = * nodeMesh2;
}

void ReadCgnsZoneAttribute( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    ONEFLOW::ReadCgnsZoneType( myZone, cgnsZoneIn );

    ONEFLOW::ReadCgnsZoneNameAndGeneralizedDimension( myZone, cgnsZoneIn );

    ONEFLOW::SetDimension( myZone, cgnsZoneIn );
}

void ReadCgnsGrid( CgnsZone * myZone, CgnsZone * cgnsZoneIn )
{
    ONEFLOW::ReadCgnsZoneAttribute( myZone, cgnsZoneIn );

    ONEFLOW::ReadElementConnectivities( myZone, cgnsZoneIn );

    ONEFLOW::ReadCgnsGridCoordinates( myZone, cgnsZoneIn );

    myZone->ConvertToInnerDataStandard();
}

void DumpCgnsZoneType( CgnsZone * myZone, Grid * grid )
{
    if ( IsUnsGrid( grid->type ) )
    {
        myZone->cgnsZoneType = CGNS_ENUMV( Unstructured );
    }
    else
    {
        myZone->cgnsZoneType = CGNS_ENUMV( Structured );
    }

    cout << "   The Zone Type is " << GetCgnsZoneTypeName( myZone->cgnsZoneType ) << " Zone" << "\n";
}

void FillISize( CgInt *isize, int ni, int nj, int nk, int dimension )
{
    int j = 0;
    // vertex size
    isize[ j ++ ] = ni;
    isize[ j ++ ] = nj;
    if ( dimension == THREE_D )
    {
        isize[ j ++ ] = nk;
    }
    // cell size
    isize[ j ++ ] = ni - 1;
    isize[ j ++ ] = nj - 1;
    if ( dimension == THREE_D )
    {
        //isize[ j ++ ] = MAX( nk - 1, 1 );
        isize[ j ++ ] = nk - 1;
    }
    // boundary vertex size (always zero for structured grids)
    isize[ j ++ ] = 0;
    isize[ j ++ ] = 0;
    if ( dimension == THREE_D )
    {
        isize[ j ++ ] = 0;
    }
}

void FillISize( CgnsZone * myZone, Grid * gridIn )
{
    StrGrid * grid = ONEFLOW::StrGridCast( gridIn );
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;
    ONEFLOW::FillISize( myZone->isize, ni, nj, nk, THREE_D );
}

void DumpCgnsZoneNameAndGeneralizedDimension( CgnsZone * myZone, Grid * gridIn )
{
    ONEFLOW::FillISize( myZone, gridIn );

    myZone->zoneName = gridIn->name;
    myZone->zId = -1;
    cout << " cell dim = " << myZone->cgnsBase->celldim << " physics dim = " << myZone->cgnsBase->phydim << "\n";
    //create zone
    cg_zone_write( myZone->cgnsBase->fileId, myZone->cgnsBase->baseId, myZone->zoneName.c_str(), myZone->isize, myZone->cgnsZoneType, &myZone->zId );
    cout << " Zone Id = " << myZone->zId << "\n";

    cout << "   CGNS Zone Name = " << myZone->zoneName << "\n";
}

void DumpCgnsZoneAttribute( CgnsZone * myZone, Grid * grid )
{
    ONEFLOW::DumpCgnsZoneType( myZone, grid );

    ONEFLOW::DumpCgnsZoneNameAndGeneralizedDimension( myZone, grid );
}

void DumpCgnsGridBoundary( CgnsZone * myZone, Grid * grid )
{
    myZone->cgnsZbc->DumpCgnsGridBoundary( grid );
}

void DumpCgnsGridCoordinates( CgnsZone * myZone, Grid * grid )
{
    // write grid coordinates (user must use SIDS-standard names here)
    int index_x = -1;
    int index_y = -2;
    int index_z = -3;
    cg_coord_write( myZone->cgnsBase->fileId, myZone->cgnsBase->baseId, myZone->zId, RealDouble, "CoordinateX", &grid->nodeMesh->xN[0], &index_x );
    cg_coord_write( myZone->cgnsBase->fileId, myZone->cgnsBase->baseId, myZone->zId, RealDouble, "CoordinateY", &grid->nodeMesh->yN[0], &index_y );
    cg_coord_write( myZone->cgnsBase->fileId, myZone->cgnsBase->baseId, myZone->zId, RealDouble, "CoordinateZ", &grid->nodeMesh->zN[0], &index_z );
    cout << " index_x = " << index_x << "\n";
    cout << " index_y = " << index_y << "\n";
    cout << " index_z = " << index_z << "\n";
}

void DumpCgnsZone( CgnsZone * myZone, Grid * grid )
{
    ONEFLOW::DumpCgnsZoneAttribute( myZone, grid );

    ONEFLOW::DumpCgnsGridBoundary( myZone, grid );

    ONEFLOW::DumpCgnsGridCoordinates( myZone, grid );
}

void PrepareCgnsZone( CgnsZone * myZone, Grid * grid )
{
    Grids grids;
    grids.push_back( grid );
    myZone->cgnsZoneType = CGNS_ENUMV( Unstructured );
    ONEFLOW::PrepareCgnsZoneSub( grids, myZone );
}



#endif
EndNameSpace