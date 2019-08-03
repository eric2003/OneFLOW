/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "CgnsBcRegion.h"
#include "CgnsBcInterface.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "Boundary.h"
#include "BasicIO.h"
#include "Stop.h"
#include "Dimension.h"
#include "HXMath.h"
#include "FaceSolver.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

CgnsBcRegion::CgnsBcRegion( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->bcInterface = 0;
}

CgnsBcRegion::~CgnsBcRegion()
{
    delete this->bcInterface;
}

void CgnsBcRegion::ConvertToInnerDataStandard()
{
    //对于各种情况均成立
    for ( int eId = 0; eId < this->nElements; ++ eId )
    {
        this->connList[ eId ] -= 1;
    }

    if ( bcInterface )
    {
        bcInterface->ConvertToInnerDataStandard();
    }

}

int CgnsBcRegion::ComputeBase()
{
    for ( int eId = 0; eId < this->nElements; ++ eId )
    {
        if ( this->connList[ eId ] < this->cgnsZone->nCell )
        {
            return 0;
        }
    }
    return 1;
}

void CgnsBcRegion::ShiftBcRegion()
{
    if ( this->gridLocation != Vertex )
    {
        for ( int eId = 0; eId < this->nElements; ++ eId )
        {
            this->connList[ eId ] += this->cgnsZone->nCell; //此处增加了偏移量，则相应的单元编号也应增加偏移量
        }
    }
}

void CgnsBcRegion::ProcessVertexBc( IntSet & bcVertex )
{
    if ( this->pointSetType == ElementRange || this->pointSetType == PointRange )
    {
        for ( int iBcPoint = this->connList[ 0 ]; iBcPoint <= this->connList[ 1 ]; ++ iBcPoint )
        {
            bcVertex.insert( this->cgnsZone->l2g[ iBcPoint ] );
        }
    }
    else
    {
        for ( int iBcPoint = 0; iBcPoint < this->nElements; ++ iBcPoint )
        {
            bcVertex.insert( this->cgnsZone->l2g[ this->connList[ iBcPoint ] ] );
        }
    }
}

void CgnsBcRegion::ProcessFaceBc( IntSet & bcVertex )
{
    if ( this->pointSetType == ElementRange ||
         this->pointSetType == PointRange )
    {
        cout << " nBcElement = " << this->connList[ 1 ] - this->connList[ 0 ] + 1 << endl;
        for ( int eId = this->connList[ 0 ]; eId <= this->connList[ 1 ]; ++ eId )
        {
            CgIntField fNodeId;
            cgnsZone->GetElementNodeId( eId, fNodeId );
            for ( int iNode = 0; iNode < fNodeId.size(); ++ iNode )
            {
                bcVertex.insert( this->cgnsZone->l2g[ fNodeId[ iNode ] ] );
            }
        }
    }
    else // if (  this->pointSetType == ElementList ||  this->pointSetType == PointList )
    {
        for ( int eId = 0; eId < this->nElements; ++ eId )
        {
            CgIntField fNodeId;
            cgnsZone->GetElementNodeId( this->connList[ eId ], fNodeId );
            for ( int iNode = 0; iNode < fNodeId.size(); ++ iNode )
            {
                bcVertex.insert( this->cgnsZone->l2g[ fNodeId[ iNode ] ] );
            }
        }
    }
}


void CgnsBcRegion::ScanBcFace( FaceSolver * face_solver )
{
    IntSet bcVertex;
    if ( this->gridLocation == Vertex )
    {
        this->ProcessVertexBc( bcVertex );
        //this->ProcessFaceBc( bcVertex );
    }
    else
    {
        this->ProcessFaceBc( bcVertex );
    }
    //face_solver->ScanBcFace( bcVertex, this->bcType, this->nameId );
    face_solver->ScanBcFaceDetail( bcVertex, this->bcType, this->nameId );
}

void CgnsBcRegion::ReadCgnsOrdinaryBcRegion()
{
    this->ReadCgnsOrdinaryBcRegionInfo();

    this->ReadCgnsOrdinaryBcRegionGridLocation();

    this->CreateCgnsBcConn();

    this->ReadCgnsBcConn();

    this->PrintCgnsBcConn();
}


void CgnsBcRegion::ReadCgnsOrdinaryBcRegionInfo()
{
    // Read the info for this boundary condition.
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    CgnsTraits::char33 bcRegionName;

    this->gridConnType = GridConnectivityTypeNull;

    cg_boco_id( fileId, baseId, zId, this->id, & this->bc_double_id );

    cg_boco_info( fileId, baseId, zId, this->id,
                  bcRegionName, & this->bcType, & this->pointSetType, & this->nElements,
                  normalIndex,  & normalListSize, & this->normalDataType, & this->nDataSets );

    this->name = bcRegionName;

    cgnsZone->cgnsBase->SetFamilyBc( bcType, bcRegionName );

    cout << "   CGNS Boundary Name             = " << bcRegionName << "\n";
    cout << "   CGNS Boundary Condition Name   = " << GetCgnsBcName( this->bcType ) << "\n";

    //if ( this->name == "outflow" && this->bcType == BCTypeNull )
    //{
    //    this->bcType = BCOutflow;
    //}
}

void CgnsBcRegion::ReadCgnsOrdinaryBcRegionGridLocation()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    cg_goto( fileId, baseId, "Zone_t", zId, "ZoneBC_t", 1, "BC_t", this->id, "end" );

    GridLocation_t bcGridLocation;
    cg_gridlocation_read( & bcGridLocation );
    this->gridLocation = bcGridLocation;

    cout << "   CGNS Grid Location Name        = " << GetCgnsGridLocationName( bcGridLocation ) << "\n";
    if ( bcGridLocation == CGNS_ENUMV( FaceCenter ) )
    {
        //cout << "   GridLocation == CGNS_ENUMV( FaceCenter ) means BC data refers to elements, not nodes\n";
    }
    else if ( bcGridLocation == CGNS_ENUMV( Vertex ) )
    {
        //cout << "   GridLocation == Vertex means BC data refers to nodes, not elements\n";

        //cout << " Now Change Vertex to FaceCenter for incompatible reasons\n";
        //iGridLocation = CGNS_ENUMV( FaceCenter );
    }
    //cout << "\n";
}

void CgnsBcRegion::CreateCgnsBcConn()
{
    //cout << "   CGNS Zone Type Name            = " << GetCgnsZoneTypeName( cgnsZone->cgnsZoneType ) << "\n";

    if ( cgnsZone->cgnsZoneType == Unstructured )
    {
        this->connList.resize( nElements );
    }
    else
    {
        if ( this->nElements == 2 )
        {
            this->nElements = 6;
        }
        this->connList.resize( nElements );
    }
}

void CgnsBcRegion::ReadCgnsBcConn()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    cout << "   CGNS PointSet Type Name        = " << GetCgnsPointSetName( this->pointSetType ) << "\n";

    int cgnsNormalList;

    // Read the element ID’s.
    cg_boco_read( fileId, baseId, zId, this->id, & connList[ 0 ], & cgnsNormalList );
    int kkk = 1;
}

void CgnsBcRegion::ProcessCgns1to1BcRegion( int bcId )
{
    this->id = bcId;
    this->bcInterface = new CgnsBcInterface( this );
    this->bcInterface->ReadCgnsBcConnInfo();
    this->bcInterface->ReadCgnsBcConnData();

    this->nElements = this->bcInterface->nConnPoints;
    this->bcType = CGNS_ENUMV( BCTypeNull );
    //this->bcType = CGNS_ENUMV( BCSymmetryPlane );

    cout << "   CGNS Boundary Name             = " << this->name << "\n";
    cout << "   CGNS Boundary Condition Name   = " << GetCgnsBcName( CGNS_ENUMV( BCTypeNull ) ) << "\n";

    this->CreateCgnsBcConn();

    cout << "   CGNS PointSet Type Name        = " << GetCgnsPointSetName( this->pointSetType ) << "\n";

    // Read the element ID’s.
    for ( int iBCElement = 0; iBCElement < this->nElements; ++ iBCElement )
    {
        this->connList[ iBCElement ] = this->bcInterface->connPoint[ iBCElement ];
    }

    this->PrintCgnsBcConn();
}

void CgnsBcRegion::ReadCgns1to1BoundaryRegion( int i1to1 )
{
    this->id = i1to1;
    this->bcInterface = new CgnsBcInterface( this );
    this->bcInterface->ReadCgnsBc1To1();

    this->nElements = this->bcInterface->nConnPoints;
    this->bcType       = CGNS_ENUMV( BCTypeNull );
    this->pointSetType = CGNS_ENUMV( PointRange );
    this->gridLocation = CGNS_ENUMV( FaceCenter );
    this->CreateCgnsBcConn();

    // Read the element ID’s.
    for ( int iBCElement = 0; iBCElement < this->nElements; ++ iBCElement )
    {
        this->connList[ iBCElement ] = this->bcInterface->connPoint[ iBCElement ];
    }

    this->PrintCgnsBcConn();
}

void CgnsBcRegion::PrintCgnsBcConn()
{
    if ( cgnsZone->cgnsZoneType == Unstructured )
    {
        if ( this->gridLocation == CGNS_ENUMV( Vertex ) )
        {
            cout << "   CGNS Boundary Point's Number   = ";
        }
        else
        {
            cout << "   CGNS Boundary Element's Number = ";
        }

        if ( this->pointSetType == CGNS_ENUMV( ElementRange ) ||
             this->pointSetType == CGNS_ENUMV( PointRange   ) )
        {
            cout << this->nElements << "( " << connList[ 1 ] - connList[ 0 ] + 1 << " )" << "\n";
        }
        else
        {
            cout << this->nElements << "\n";
        }

        if ( this->nElements == 2 )
        {
            cout << "   " << connList[ 0 ] << " " << connList[ 1 ] << "\n";
        }
        else
        {
            cout << "   ";
            CgInt num = 5;
            for ( int i = 0; i < MIN(num, this->nElements ); ++ i )
            {
                cout << connList[ i ] << " ";
            }
            cout << "\n";
        }

    }
    else
    {
        int celldim = cgnsZone->cgnsBase->celldim;

        cout << "   The Boundary Range is :\n";
        StringField rangeTitle;
        rangeTitle.push_back( "   I Direction " );
        rangeTitle.push_back( "   J Direction " );
        rangeTitle.push_back( "   K Direction " );

        IntField ijkMin( 3 ), ijkMax( 3 );

        this->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );

        for ( int iDimension = 0; iDimension < celldim; ++ iDimension )
        {
            cout << rangeTitle[ iDimension ] << setw( 10 ) << ijkMin[ iDimension ] << setw( 10 ) <<  ijkMax[ iDimension ] << "\n";
        }
        cout << "\n";
    }
}

void CgnsBcRegion::ReconstructStrRegion( IntField & ijkMin, IntField & ijkMax )
{
    this->bcType       = CGNS_ENUMV( BCInflow );
    this->pointSetType = CGNS_ENUMV( PointRange );
    this->gridLocation = CGNS_ENUMV( FaceCenter );
    this->nElements = 6;
    this->connList.resize( this->nElements );

    string bcName = GetCgnsBcName( this->bcType );
    this->name = AddString( this->cgnsZone->zoneName, bcName, this->id );

    int celldim = cgnsZone->cgnsBase->celldim;

    if ( celldim == TWO_D )
    {
        this->connList[ 0 ] = ijkMin[ 0 ];
        this->connList[ 1 ] = ijkMin[ 1 ];

        this->connList[ 2 ] = ijkMax[ 0 ];
        this->connList[ 3 ] = ijkMax[ 1 ];
    }
    else
    {
        this->connList[ 0 ] = ijkMin[ 0 ];
        this->connList[ 1 ] = ijkMin[ 1 ];
        this->connList[ 2 ] = ijkMin[ 2 ];

        this->connList[ 3 ] = ijkMax[ 0 ];
        this->connList[ 4 ] = ijkMax[ 1 ];
        this->connList[ 5 ] = ijkMax[ 2 ];
    }
 }

void CgnsBcRegion::ExtractIJKRegionFromBcConn( IntField & ijkMin, IntField & ijkMax )
{
    this->ExtractIJKRegionFromBcConn( ijkMin, ijkMax, this->connList );
}

void CgnsBcRegion::ExtractIJKRegionFromBcConn( IntField & ijkMin, IntField & ijkMax, CgIntField& bcConn )
{
    int imin, imax, jmin, jmax, kmin, kmax;
    int celldim = cgnsZone->cgnsBase->celldim;
    if ( celldim == TWO_D )
    {
        imin = bcConn[ 0 ];
        jmin = bcConn[ 1 ];
        kmin = 1;
        imax = bcConn[ 2 ];
        jmax = bcConn[ 3 ];
        kmax = 1;
    }
    else
    {
        imin = bcConn[ 0 ];
        jmin = bcConn[ 1 ];
        kmin = bcConn[ 2 ];
        imax = bcConn[ 3 ];
        jmax = bcConn[ 4 ];
        kmax = bcConn[ 5 ];
    }

    ijkMin[ 0 ] = MIN( ABS( imin ), ABS( imax ) );
    ijkMin[ 1 ] = MIN( ABS( jmin ), ABS( jmax ) );
    ijkMin[ 2 ] = MIN( ABS( kmin ), ABS( kmax ) );

    ijkMax[ 0 ] = MAX( ABS( imin ), ABS( imax ) );
    ijkMax[ 1 ] = MAX( ABS( jmin ), ABS( jmax ) );
    ijkMax[ 2 ] = MAX( ABS( kmin ), ABS( kmax ) );
}

void CgnsBcRegion::CopyStrBcRegion( CgnsBcRegion * strBcRegion, CgInt & startId )
{
    this->name = strBcRegion->name;
    this->nElements = 2;
    this->bcType = strBcRegion->bcType;
    this->pointSetType = CGNS_ENUMV( ElementRange );
    this->gridLocation = CGNS_ENUMV( CellCenter   );

    this->CreateCgnsBcConn();

    this->ReadCgnsBcConn( strBcRegion, startId );
}

void CgnsBcRegion::ReadCgnsBcConn( CgnsBcRegion * strBcRegion, CgInt& startId )
{
    CgInt actualNumberOfBoundaryElement = strBcRegion->GetActualNumberOfBoundaryElements();
    this->connList[ 0 ] = startId;
    this->connList[ 1 ] = actualNumberOfBoundaryElement - 1 + startId;
    startId += actualNumberOfBoundaryElement;
}

CgInt CgnsBcRegion::GetActualNumberOfBoundaryElements()
{
    if ( cgnsZone->cgnsZoneType == Unstructured )
    {
        if ( this->gridLocation == CGNS_ENUMV( Vertex ) )
        {
            cout << "   Can't Get Element Number For cgnsBoundaryGridLocation == Vertex Case \n";
            return INVALID_INDEX;
        }

        if ( this->pointSetType == CGNS_ENUMV( ElementRange ) ||
             this->pointSetType == CGNS_ENUMV( PointRange   ) )
        {
            return this->connList[ 1 ] - this->connList[ 0 ] + 1;
        }
        else
        {
            return this->nElements;
        }
    }
    else
    {
        int celldim = cgnsZone->cgnsBase->celldim;

        IntField ijkMin( 3 ), ijkMax( 3 );
        this->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );

        int actualNumberOfBoundaryElements = 1;

        for ( int idim = 0; idim < celldim; ++ idim )
        {
            actualNumberOfBoundaryElements *= MAX( ijkMax[ idim ] - ijkMin[ idim ], 1 );
        }
        return actualNumberOfBoundaryElements;
    }
}

void SetBcConn( CgnsZone * cgnsZone, IntField & ijkMin, IntField & ijkMax, CgIntField& conn, int & pos, int & nElem )
{
    int ni = static_cast<int> (cgnsZone->GetNI());
    int nj = static_cast<int> (cgnsZone->GetNJ());
    int nk = static_cast<int> (cgnsZone->GetNK());

    int ist, jst, kst, ied, jed, ked;

    ist = ijkMin[ 0 ];
    jst = ijkMin[ 1 ];
    kst = ijkMin[ 2 ];

    ied = ijkMax[ 0 ];
    jed = ijkMax[ 1 ];
    ked = ijkMax[ 2 ];

    cout << " ist, ied, jst, jed, kst, ked = " << ist << " " << ied << " " << jst << " " << jed << " " << kst << " " << ked << "\n";

    int celldim = cgnsZone->cgnsBase->celldim;
    int numpt = 4;
    if ( celldim == TWO_D ) numpt = 2;
    if ( celldim == ONE_D ) numpt = 1;

    int index1, index2, index3, index4;

    if ( celldim == ONE_D )
    {
        int i = ist;
        int j = jst;
        int k = kst;
        EncodeIJK( index1, i, j, k, ni, nj, nk );
        conn[ pos + 0 ] = index1 + 1;
        pos += numpt;
        return;
    }

    if ( ist == ied )
    {
        int i = ist;
        if ( celldim == THREE_D )
        {
            for ( int k = kst; k <= ked - 1; ++ k )
            {
                for ( int j = jst; j <= jed - 1; ++ j )
                {
                    if ( i == 1 )
                    {
                        EncodeIJK( index1,  i, j    , k    ,  ni, nj, nk );
                        EncodeIJK( index2,  i, j    , k + 1,  ni, nj, nk );
                        EncodeIJK( index3,  i, j + 1, k + 1,  ni, nj, nk );
                        EncodeIJK( index4,  i, j + 1, k    ,  ni, nj, nk );
                    }
                    else
                    {
                        EncodeIJK( index1,  i, j    , k    ,  ni, nj, nk );
                        EncodeIJK( index2,  i, j + 1, k    ,  ni, nj, nk );
                        EncodeIJK( index3,  i, j + 1, k + 1,  ni, nj, nk );
                        EncodeIJK( index4,  i, j    , k + 1,  ni, nj, nk );
                    }

                    conn[ pos + 0 ] = index1 + 1;
                    conn[ pos + 1 ] = index2 + 1;
                    conn[ pos + 2 ] = index3 + 1;
                    conn[ pos + 3 ] = index4 + 1;

                    pos += numpt;
                    ++ nElem;
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
                    EncodeIJK( index1,  i, j + 1, k    ,  ni, nj, nk );
                    EncodeIJK( index2,  i, j    , k    ,  ni, nj, nk );
                }
                else
                {
                    EncodeIJK( index1,  i, j    , k    ,  ni, nj, nk );
                    EncodeIJK( index2,  i, j + 1, k    ,  ni, nj, nk );
                }
                conn[ pos + 0 ] = index1 + 1;
                conn[ pos + 1 ] = index2 + 1;
                pos += numpt;
            }
        }
        return;
    }

    if ( jst == jed )
    {
        int j = jst;
        if ( celldim == THREE_D )
        {
            for ( int k = kst; k <= ked - 1; ++ k )
            {
                for ( int i = ist; i <= ied - 1; ++ i )
                {
                    if ( j == 1 )
                    {
                        EncodeIJK( index1, i    , j   , k    ,  ni, nj, nk );
                        EncodeIJK( index2, i + 1, j   , k    ,  ni, nj, nk );
                        EncodeIJK( index3, i + 1, j   , k + 1,  ni, nj, nk );
                        EncodeIJK( index4, i    , j   , k + 1,  ni, nj, nk );
                    }   
                    else
                    {
                        EncodeIJK( index1, i    , j    , k    ,  ni, nj, nk );
                        EncodeIJK( index2, i    , j    , k + 1,  ni, nj, nk );
                        EncodeIJK( index3, i + 1, j    , k + 1,  ni, nj, nk );
                        EncodeIJK( index4, i + 1, j    , k    ,  ni, nj, nk );
                    }
                    conn[ pos + 0 ] = index1 + 1;
                    conn[ pos + 1 ] = index2 + 1;
                    conn[ pos + 2 ] = index3 + 1;
                    conn[ pos + 3 ] = index4 + 1;

                    pos += numpt;
                    ++ nElem;
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
                    EncodeIJK( index1, i    , j   , k    ,  ni, nj, nk );
                    EncodeIJK( index2, i + 1, j   , k    ,  ni, nj, nk );
                }   
                else
                {
                    EncodeIJK( index1, i + 1, j    , k   ,  ni, nj, nk );
                    EncodeIJK( index2, i    , j    , k   ,  ni, nj, nk );
                }

                conn[ pos + 0 ] = index1 + 1;
                conn[ pos + 1 ] = index2 + 1;
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
                    EncodeIJK( index1, i    , j    , k    ,  ni, nj, nk );
                    EncodeIJK( index2, i    , j + 1, k    ,  ni, nj, nk );
                    EncodeIJK( index3, i + 1, j + 1, k    ,  ni, nj, nk );
                    EncodeIJK( index4, i + 1, j    , k    ,  ni, nj, nk );
                }   
                else
                {
                    EncodeIJK( index1, i    , j    , k    ,  ni, nj, nk );
                    EncodeIJK( index2, i + 1, j    , k    ,  ni, nj, nk );
                    EncodeIJK( index3, i + 1, j + 1, k    ,  ni, nj, nk );
                    EncodeIJK( index4, i    , j + 1, k    ,  ni, nj, nk );
                }

                conn[ pos + 0 ] = index1 + 1;
                conn[ pos + 1 ] = index2 + 1;
                conn[ pos + 2 ] = index3 + 1;
                conn[ pos + 3 ] = index4 + 1;

                pos += numpt;
                ++ nElem;
            }
        }

        return;
    }

    Stop( " error : ist != ied, jst != jed, kst != ked \n" );
}

#endif

EndNameSpace