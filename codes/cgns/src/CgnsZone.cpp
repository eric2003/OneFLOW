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

CgnsZone::CgnsZone( CgnsBase * cgnsBase )
{
    this->cgnsBase = cgnsBase;
    this->cgnsZsection = 0;
    this->cgnsZbc = 0;
    this->volBcType = -1;
    this->cgnsCoor = 0;
}

CgnsZone::~CgnsZone()
{
    delete this->cgnsZsection;
    delete this->cgnsZbc;
    delete this->cgnsCoor;
}

void CgnsZone::SetVolBcType( int volBcType )
{
    this->volBcType = volBcType;
}

int CgnsZone::GetVolBcType()
{
    return this->volBcType;
}

void CgnsZone::Create()
{
    this->cgnsZsection = new CgnsZsection( this );
    this->cgnsZbc = new CgnsZbc( this );
    this->cgnsCoor = new CgnsCoor( this );
}

void CgnsZone::InitElement( GridElem * ge )
{
    this->ConstructCgnsGridPoints( ge->point_factory );
    this->SetElementTypeAndNode( ge->elem_feature );
}

void CgnsZone::SetPeriodicBc()
{
    this->cgnsZbc->SetPeriodicBc();
}

void CgnsZone::SetElementTypeAndNode( ElemFeature * elem_feature )
{
    size_t nSection = this->cgnsZsection->nSection;
    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->cgnsZsection->GetCgnsSection( iSection );
        cgnsSection->SetElementTypeAndNode( elem_feature );
    }
    cout << "\n";
    cout << " iZone = " << this->zId << " nCell = " << this->cgnsCoor->GetNCell() << "\n";
    cout << " elem_feature->eType->size = " << elem_feature->eType->size() << endl;
}

bool CgnsZone::ExistSection( const string & sectionName )
{
    return this->cgnsZsection->ExistSection( sectionName );
}

void CgnsZone::InitLgMapping()
{
    int nNode = this->cgnsCoor->GetNNode();
    this->l2g.resize( nNode );

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        this->l2g[ iNode ] = iNode;
    }
}

void CgnsZone::ConvertToInnerDataStandard()
{
    if ( this->cgnsZoneType == CGNS_ENUMV( Structured ) )
    {
        //this->SetStructuredSectionInformation();
        return;
    }

    this->cgnsZsection->ConvertToInnerDataStandard();
    this->cgnsZbc->ConvertToInnerDataStandard();

}

void CgnsZone::ConstructCgnsGridPoints( PointFactory * point_factory )
{
    NodeMesh * nodeMesh = this->cgnsCoor->GetNodeMesh();
    RealField & x = nodeMesh->xN;
    RealField & y = nodeMesh->yN;
    RealField & z = nodeMesh->zN;

    this->InitLgMapping();

    size_t nNode = nodeMesh->GetNumberOfNodes();

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        int pid = point_factory->AddPoint( x[ iNode ], y[ iNode ], z[ iNode ] );
        
        this->l2g[ iNode ] = pid; //For the merging of partitioned multiple blocks with single block
    }
}

void CgnsZone::ScanBcFace( FaceSolver * face_solver )
{
    this->cgnsZbc->ScanBcFace( face_solver );
}

void CgnsZone::GetElementNodeId( CgInt eId, CgIntField & eNodeId )
{
    CgnsSection * cgnsSection = this->cgnsZsection->GetSectionByEid( eId );
    cgnsSection->GetElementNodeId( eId - cgnsSection->startId, eNodeId );
}

void CgnsZone::FillISize( Grid * gridIn )
{
    StrGrid * grid = ONEFLOW::StrGridCast( gridIn );
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;
    this->FillISize( ni, nj, nk, THREE_D );
}

void CgnsZone::FillISize( int ni, int nj, int nk, int dimension )
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

void CgnsZone::DumpCgnsZone( Grid * grid )
{
    this->DumpCgnsZoneAttribute( grid );

    this->DumpCgnsGridBoundary( grid );

    this->DumpCgnsGridCoordinates( grid );
}

void CgnsZone::ReadCgnsGrid()
{
    this->ReadCgnsZoneAttribute();

    this->ReadElementConnectivities();

    this->ReadCgnsGridBoundary();

    this->ReadCgnsGridCoordinates();

    this->ConvertToInnerDataStandard();
}

void CgnsZone::ReadCgnsGrid( CgnsZone * cgnsZoneIn )
{
    this->ReadCgnsZoneAttribute( cgnsZoneIn );

    this->ReadElementConnectivities( cgnsZoneIn );

    this->ReadCgnsGridCoordinates( cgnsZoneIn );

    this->ConvertToInnerDataStandard();
}

void CgnsZone::ReadCgnsZoneAttribute()
{
    this->ReadCgnsZoneType();

    this->ReadCgnsZoneNameAndGeneralizedDimension();

    this->SetDimension();
}

void CgnsZone::DumpCgnsZoneAttribute( Grid * grid )
{
    this->DumpCgnsZoneType( grid );

    this->DumpCgnsZoneNameAndGeneralizedDimension( grid );
}

void CgnsZone::ReadCgnsZoneAttribute( CgnsZone * cgnsZoneIn )
{
    this->ReadCgnsZoneType( cgnsZoneIn );

    this->ReadCgnsZoneNameAndGeneralizedDimension( cgnsZoneIn );

    this->SetDimension( cgnsZoneIn );
}

void CgnsZone::ReadCgnsZoneType()
{
    //Check the zone type
    cg_zone_type( cgnsBase->fileId, cgnsBase->baseId, this->zId, & cgnsZoneType );

    cout << "   The Zone Type is " << GetCgnsZoneTypeName( cgnsZoneType ) << " Zone" << "\n";
}

void CgnsZone::DumpCgnsZoneType( Grid * grid )
{
    if ( IsUnsGrid( grid->type ) )
    {
        this->cgnsZoneType = CGNS_ENUMV( Unstructured );
    }
    else
    {
        this->cgnsZoneType = CGNS_ENUMV( Structured );
    }

    cout << "   The Zone Type is " << GetCgnsZoneTypeName( cgnsZoneType ) << " Zone" << "\n";
}

void CgnsZone::ReadCgnsZoneType( CgnsZone * cgnsZoneIn )
{
    this->cgnsZoneType = CGNS_ENUMV( Unstructured );

    cout << "   The Zone Type is " << GetCgnsZoneTypeName( cgnsZoneType ) << " Zone" << "\n";
}

void CgnsZone::ReadCgnsZoneNameAndGeneralizedDimension()
{
    CgnsTraits::char33 cgnsZoneName;

    //Determine the number of vertices and cellVolume elements in this zone
    cg_zone_read( cgnsBase->fileId, cgnsBase->baseId, this->zId, cgnsZoneName, this->isize );

    this->zoneName = cgnsZoneName;

    cout << "   CGNS Zone Name = " << cgnsZoneName << "\n";
}

void CgnsZone::DumpCgnsZoneNameAndGeneralizedDimension( Grid * gridIn )
{
    this->FillISize( gridIn );

    this->zoneName = gridIn->name;
    this->zId = -1;
    cout << " cell dim = " << this->cgnsBase->celldim << " physics dim = " << this->cgnsBase->phydim << "\n";
    //create zone
    cg_zone_write(cgnsBase->fileId, cgnsBase->baseId, this->zoneName.c_str(), this->isize, this->cgnsZoneType, &this->zId );
    cout << " Zone Id = " << this->zId << "\n";

    cout << "   CGNS Zone Name = " << zoneName << "\n";
}

void CgnsZone::ReadCgnsZoneNameAndGeneralizedDimension( CgnsZone * cgnsZoneIn )
{
    this->zoneName = cgnsZoneIn->zoneName;
}

void CgnsZone::SetDimension()
{
    this->cgnsCoor->SetDimension();
    //int rind[ 6 ];
    //int result = cg_rind_read( rind );
    //if ( result != CG_OK )
    //{
    //    for ( int i = 0; i < 6; ++ i )
    //    {
    //        rind[i] = 0;
    //    }
    //}

    //if ( this->cgnsZoneType == CGNS_ENUMV( Structured ) )
    //{
    //    // lower range index
    //    irmin[ 0 ] = 1;
    //    irmin[ 1 ] = 1;
    //    irmin[ 2 ] = 1;

    //    // upper range index of vertices
    //    irmax[ 0 ] = 1;
    //    irmax[ 1 ] = 1;
    //    irmax[ 2 ] = 1;

    //    cellSize[ 0 ] = 1;
    //    cellSize[ 1 ] = 1;
    //    cellSize[ 2 ] = 1;

    //    // upper range index of vertices
    //    // vertex size
    //    int j = 0;
    //    irmax[ 0 ] = isize[ j ++ ];
    //    irmax[ 1 ] = isize[ j ++ ];
    //    if ( this->cgnsBase->celldim == THREE_D )
    //    {
    //        irmax[ 2 ] = isize[ j ++ ];
    //    }
    //    // cell size
    //    cellSize[ 0 ] = isize[ j ++ ];
    //    cellSize[ 1 ] = isize[ j ++ ];
    //    if ( this->cgnsBase->celldim == THREE_D )
    //    {
    //        cellSize[ 2 ] = isize[ j ++ ];
    //    }
    //    cout << "   The Dimension Of Grid is : \n";
    //    cout << "   I Direction " << setw( 10 ) << irmin[ 0 ] << setw( 10 ) << irmax[ 0 ] << "\n";
    //    cout << "   J Direction " << setw( 10 ) << irmin[ 1 ] << setw( 10 ) << irmax[ 1 ] << "\n";
    //    cout << "   K Direction " << setw( 10 ) << irmin[ 2 ] << setw( 10 ) << irmax[ 2 ] << "\n";
    //    int nNode = irmax[ 0 ] * irmax[ 1 ] * irmax[ 2 ];
    //    int nCell = cellSize[ 0 ] * cellSize[ 1 ] * cellSize[ 2 ];
    //    this->cgnsCoor->SetNNode( nNode );
    //    this->cgnsCoor->SetNCell( nCell );
    //}
    //else
    //{
    //    irmin[ 0 ] = 1;
    //    irmin[ 1 ] = 0;
    //    irmin[ 2 ] = 0;

    //    irmax[ 0 ] = isize[ 0 ];
    //    irmax[ 1 ] = 0;
    //    irmax[ 2 ] = 0;

    //    cellSize[ 0 ] = isize[ 1 ];

    //    int nNode = irmax[ 0 ];
    //    int nCell = cellSize[ 0 ];
    //    this->cgnsCoor->SetNNode( nNode );
    //    this->cgnsCoor->SetNCell( nCell );

    //}

    //cout << "   numberOfNodes = " << this->cgnsCoor->GetNNode() << " numberOfCells = " << this->cgnsCoor->GetNCell() << "\n";
}

void CgnsZone::SetDimension( CgnsZone * cgnsZoneIn )
{
    CgnsCoor * cgnsCoorIn = cgnsZoneIn->cgnsCoor;
    this->cgnsCoor->SetDimension( cgnsCoorIn );
    //isize[ 0 ] = cgnsZoneIn->cgnsCoor->GetNNode();
    //isize[ 1 ] = cgnsZoneIn->cgnsCoor->GetNCell();

    //irmin[ 0 ] = 1;
    //irmin[ 1 ] = 0;
    //irmin[ 2 ] = 0;

    //irmax[ 0 ] = isize[ 0 ];
    //irmax[ 1 ] = 0;
    //irmax[ 2 ] = 0;

    //cellSize[ 0 ] = isize[ 1 ];

    //this->cgnsCoor->SetNNode( irmax[ 0 ] );
    //this->cgnsCoor->SetNCell( cellSize[ 0 ] );

    //this->InitL2g();

    //cout << "   numberOfNodes = " << this->cgnsCoor->GetNNode() << " numberOfCells = " << this->cgnsCoor->GetNCell() << "\n";
}

void CgnsZone::InitL2g()
{
    int nNode = this->cgnsCoor->GetNNode();
    l2g.resize( nNode );

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        l2g[ iNode ] = iNode;
    }
}

CgInt CgnsZone::GetNI() const { return this->cgnsCoor->irmax[0]; };
CgInt CgnsZone::GetNJ() const { return this->cgnsCoor->irmax[1]; };
CgInt CgnsZone::GetNK() const { return this->cgnsCoor->irmax[2]; };

void CgnsZone::ReadElementConnectivities()
{
    if ( this->cgnsZoneType == CGNS_ENUMV( Structured ) ) return;

    this->ReadNumberOfCgnsSections();

    this->CreateCgnsSections();

    this->ReadCgnsSections();
}

void CgnsZone::ReadElementConnectivities( CgnsZone * cgnsZoneIn )
{
    this->AllocateUnsElemConn   ( cgnsZoneIn );
    this->GenerateUnsVolElemConn( cgnsZoneIn );
    this->GenerateUnsBcElemConn ( cgnsZoneIn );
    this->SetElemPosition();
    this->GenerateUnsBcCondConn ( cgnsZoneIn );
}

void CgnsZone::AllocateUnsElemConn( CgnsZone * cgnsZoneIn )
{
    this->cgnsZsection->nSection = 2;
    this->cgnsZsection->CreateCgnsSection();

    int s1, e1, s2, e2, etype1, etype2;
    cgnsZoneIn->GetStrZonePara( s1, e1, s2, e2, etype1, etype2 );

    CgnsSection * cgnsSection1 = this->cgnsZsection->GetCgnsSection( 0 );
    CgnsSection * cgnsSection2 = this->cgnsZsection->GetCgnsSection( 1 );
    cgnsSection1->SetSectionInfo( "Section1", etype1, s1, e1 );
    cgnsSection2->SetSectionInfo( "Section2", etype2, s2, e2 );

    this->cgnsZsection->CreateConnList();
}

void CgnsZone::GetStrZonePara( int & s1, int & e1, int & s2, int & e2, int & etype1, int & etype2  )
{
   int nActualBcFace = this->cgnsZbc->GetNumberOfActualBcElements();

    s1 = 1;
    e1 = this->cgnsCoor->GetNCell();

    s2 = e1 + 1;
    e2 = e1 + nActualBcFace;

    int celldim = this->cgnsBase->celldim;

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

void CgnsZone::SetElemPosition()
{
    this->cgnsZsection->SetElemPosition();
}

void CgnsZone::GenerateUnsVolElemConn( CgnsZone * cgnsZoneIn )
{
    int ni = static_cast<int> (cgnsZoneIn->GetNI());
    int nj = static_cast<int> (cgnsZoneIn->GetNJ());
    int nk = static_cast<int> (cgnsZoneIn->GetNK());

    cout << " ni = " << ni << " nj = " << nj << " nk = " << nk << "\n";

    int iSection = 0;
    CgnsSection * cgnsSection = this->cgnsZsection->GetCgnsSection( iSection );

    Range I, J, K;
    GetRange( ni, nj, nk, 0, -1, I, J, K );

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion( I, J, K, ist, ied, jst, jed, kst, ked );

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;

    int cell_dim = this->cgnsBase->celldim;

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

                connList[ pos ++ ] = this->l2g[ index1 ] + 1;
                connList[ pos ++ ] = this->l2g[ index2 ] + 1;

                if ( cell_dim == ONE_D ) continue;

                EncodeIJK( index3,  i + il1, j + jl1, k,  ni,  nj,  nk );
                EncodeIJK( index4,  i      , j + jl1, k,  ni,  nj,  nk );

                connList[ pos ++ ] = this->l2g[ index3 ] + 1;
                connList[ pos ++ ] = this->l2g[ index4 ] + 1;

                if ( cell_dim == TWO_D ) continue;

                int index5, index6, index7, index8;
                EncodeIJK( index5,  i      , j      , k + kl1,  ni,  nj,  nk );
                EncodeIJK( index6,  i + il1, j      , k + kl1,  ni,  nj,  nk );
                EncodeIJK( index7,  i + il1, j + jl1, k + kl1,  ni,  nj,  nk );
                EncodeIJK( index8,  i      , j + jl1, k + kl1,  ni,  nj,  nk );

                connList[ pos ++ ] = this->l2g[ index5 ] + 1;
                connList[ pos ++ ] = this->l2g[ index6 ] + 1;
                connList[ pos ++ ] = this->l2g[ index7 ] + 1;
                connList[ pos ++ ] = this->l2g[ index8 ] + 1;
            }
        }
    }
}

void CgnsZone::GenerateUnsBcElemConn( CgnsZone * cgnsZoneIn )
{
    int iSection = 1;
    CgnsSection * cgnsSection = this->cgnsZsection->GetCgnsSection( iSection );

    this->cgnsZbc->CreateCgnsZbc( cgnsZoneIn->cgnsZbc );

    cout << " ConnectionList Size = " << cgnsSection->connSize << "\n";
    cgnsZoneIn->cgnsZbc->GenerateUnsBcElemConn( cgnsSection->connList );
}

void CgnsZone::GenerateUnsBcCondConn( CgnsZone * cgnsZoneIn )
{
    int iSection = 1;
    CgnsSection * cgnsSection = this->cgnsZsection->GetCgnsSection( iSection );

    CgInt startId = cgnsSection->startId;

    int nBoco = cgnsZoneIn->cgnsZbc->cgnsZbcBoco->nBoco;
    for ( int iBoco = 0; iBoco < nBoco; ++ iBoco )
    {
        CgnsBcBoco * bcRegion    = this      ->cgnsZbc->cgnsZbcBoco->GetCgnsBc( iBoco );
        CgnsBcBoco * strBcRegion = cgnsZoneIn->cgnsZbc->cgnsZbcBoco->GetCgnsBc( iBoco );
        bcRegion->CopyStrBcRegion( strBcRegion, startId );
    }
}

void CgnsZone::ReadNumberOfCgnsSections()
{
    this->cgnsZsection->ReadNumberOfCgnsSections();
}

void CgnsZone::ReadNumberOfCgnsSections( CgnsZone * cgnsZoneIn )
{
    this->cgnsZsection->nSection = cgnsZoneIn->cgnsZsection->nSection;
    cout << "   numberOfCgnsSections = " << this->cgnsZsection->nSection << "\n";
}

void CgnsZone::CreateCgnsSections()
{
    this->cgnsZsection->CreateCgnsSection();
}

void CgnsZone::ReadCgnsSections()
{
    this->cgnsZsection->ReadCgnsSections();
}

void CgnsZone::ReadCgnsGridCoordinates()
{
    cgnsCoor->ReadCgnsGridCoordinates();
}

void CgnsZone::DumpCgnsGridCoordinates( Grid * grid )
{
    // write grid coordinates (user must use SIDS-standard names here)
    int index_x = -1;
    int index_y = -2;
    int index_z = -3;
    cg_coord_write( cgnsBase->fileId, cgnsBase->baseId, this->zId, RealDouble, "CoordinateX", &grid->nodeMesh->xN[0], &index_x );
    cg_coord_write( cgnsBase->fileId, cgnsBase->baseId, this->zId, RealDouble, "CoordinateY", &grid->nodeMesh->yN[0], &index_y );
    cg_coord_write( cgnsBase->fileId, cgnsBase->baseId, this->zId, RealDouble, "CoordinateZ", &grid->nodeMesh->zN[0], &index_z );
    cout << " index_x = " << index_x << "\n";
    cout << " index_y = " << index_y << "\n";
    cout << " index_z = " << index_z << "\n";
}

void CgnsZone::ReadCgnsGridCoordinates( CgnsZone * cgnsZoneIn )
{
    NodeMesh * nodeMesh1 = this->cgnsCoor->GetNodeMesh();
    NodeMesh * nodeMesh2 = cgnsZoneIn->cgnsCoor->GetNodeMesh();

    * nodeMesh1 = * nodeMesh2;
}

void CgnsZone::ReadCgnsGridBoundary()
{
    cgnsZbc->ReadCgnsGridBoundary();
}

void CgnsZone::DumpCgnsGridBoundary( Grid * grid )
{
    cgnsZbc->DumpCgnsGridBoundary( grid );
}

void CgnsZone::ProcessPeriodicBc()
{
    ;
}

void CgnsZone::PrepareCgnsZone( Grid * grid )
{
    Grids grids;
    grids.push_back( grid );
    this->cgnsZoneType = CGNS_ENUMV( Unstructured );
    ONEFLOW::PrepareCgnsZone( grids, this );
}


#endif
EndNameSpace