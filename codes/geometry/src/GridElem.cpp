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

#include "GridElem.h"
#include "GridTypes.h"
#include "CgnsZone.h"
#include "CgnsZbase.h"
#include "DataBase.h"
#include "HXCgns.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "CellTopo.h"
#include "CellMesh.h"
#include "FaceMesh.h"
#include "ElemFeature.h"
#include "FaceTopo.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "PointManager.h"
#include "NodeMesh.h"
#include "GridState.h"
#include "BgGrid.h"
#include "CgnsZsection.h"
#include "CgnsSection.h"
#include "Fatal.h"
#include <iostream>
#include <iomanip>
#include <utility>


BeginNameSpace( ONEFLOW )

int OneFlow2CgnsZoneType( int zoneType )
{
    if ( zoneType == UMESH )
    {
        return CGNS_ENUMV( Unstructured );
    }
    else
    {
        return CGNS_ENUMV( Structured );
    }
}

int Cgns2OneFlowZoneType( int zoneType )
{
    if ( zoneType == CGNS_ENUMV( Unstructured ) )
    {
        return UMESH;
    }
    else
    {
        return SMESH;
    }
}

GridElem::GridElem( HXVector< CgnsZone * > zoneViews )
    : zoneViews( std::move( zoneViews ) ),
      minLen( LARGE ),
      maxLen( -LARGE )
{
}

GridElem::~GridElem() = default;

CgnsZone * GridElem::GetCgnsZone( int iZone )
{
    return this->zoneViews[ iZone ];
}

const CgnsZone * GridElem::GetCgnsZone( int iZone ) const
{
    return this->zoneViews[ iZone ];
}

int GridElem::GetNZones() const
{
    return this->zoneViews.size();
}

bool GridElem::HasPolygonSection() const
{
    if ( this->GetNZones() == 0 )
    {
        Fatal( "GridElem requires at least one CGNS zone." );
    }

    const bool hasPolygon = this->GetCgnsZone( 0 )->cgnsZsection->HasPolygonSection();

    for ( int iZone = 1; iZone < this->GetNZones(); ++ iZone )
    {
        const bool zoneHasPolygon =
            this->GetCgnsZone( iZone )->cgnsZsection->HasPolygonSection();

        if ( zoneHasPolygon != hasPolygon )
        {
            Fatal( "GridElem cannot combine CGNS zones with different element-generation modes." );
        }
    }

    return hasPolygon;
}

int GridElem::GetVolBcType() const
{
    if ( this->GetNZones() == 0 )
    {
        Fatal( "GridElem requires at least one CGNS zone." );
    }

    const int volBcType = this->GetCgnsZone( 0 )->GetVolBcType();

    for ( int iZone = 1; iZone < this->GetNZones(); ++ iZone )
    {
        if ( this->GetCgnsZone( iZone )->GetVolBcType() != volBcType )
        {
            Fatal( "GridElem cannot combine CGNS zones with different volume boundary types." );
        }
    }

    return volBcType;
}

void GridElem::PrepareUnsCalcGrid()
{
    const bool flag = this->HasPolygonSection();
    if ( flag )
    {
        this->PrepareUnsCalcGridPolyhedron();
    }
    else
    {
        this->PrepareUnsCalcGridNormal();
    }
}

void GridElem::PrepareUnsCalcGridNormal()
{
    std::cout << " InitCgnsElements()\n";
    this->InitCgnsElements();
    std::cout << " ScanElements()\n";
    this->elem_feature.ScanElements( this->face_solver );
    std::cout << " ScanBcFace()\n";
    this->ScanBcFace();

    //Continue to parse
    std::cout << " ScanElements()\n";
    this->elem_feature.ScanElements( this->face_solver );
    this->GenerateCalcElement();
}

void GridElem::PrepareUnsCalcGridPolyhedron()
{
    std::cout << " ScanPolygonFace()\n";
    this->ScanPolygonFace();
    this->ScanBcFace();
    this->GenerateCalcElement();
}

void GridElem::ScanPolygonFace()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );

        cgnsZone->ConstructCgnsGridPoints( &this->point_factory );

        //Scan NGON_n PolygonFace
        int nSections = cgnsZone->cgnsZsection->nSection;
        for ( int iSection = 0; iSection < nSections; ++ iSection )
        {
            CgnsSection * cgnsSection = cgnsZone->cgnsZsection->GetCgnsSection( iSection );
            if ( cgnsSection->eType != NGON_n ) continue;
            this->face_solver.ScanPolygonFace( cgnsSection );
        }
        //Scan NFACE_n PolyhedronElement
        for ( int iSection = 0; iSection < nSections; ++ iSection )
        {
            CgnsSection * cgnsSection = cgnsZone->cgnsZsection->GetCgnsSection( iSection );
            if ( cgnsSection->eType != NFACE_n ) continue;
            this->face_solver.ScanPolyhedronElement( cgnsSection );
            this->SetPolyhedronElementType( *cgnsSection );
        }

        int nFaces = this->face_solver.faceTopo->faces.size();

    }
}

void GridElem::SetPolyhedronElementType( CgnsSection & cgnsSection )
{
    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        int e_type = cgnsSection.eTypeList[ iElem ];

        this->elem_feature.eTypes.push_back( e_type );
    }
}

void GridElem::InitCgnsElements()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        
        cgnsZone->ConstructCgnsGridPoints( &this->point_factory );
        cgnsZone->SetElementTypeAndNode( &this->elem_feature );
    }
}

void GridElem::ScanBcFace()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->ScanBcFace( this->face_solver );
    }

    this->face_solver.ScanInterfaceBc();
}

void GridElem::GenerateCalcElement()
{
    int nElement =  this->elem_feature.eTypes.size();

    FaceTopo * faceTopo = this->face_solver.faceTopo.get();

    int nFaces = this->face_solver.faceTopo->faces.size();
    int nBFaces = 0;

    //std::cout << " nFaces = " << nFaces << "\n";

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            std::cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFaces << std::endl;
        }

        int rc = ( faceTopo->rCells )[ iFace ];

        if ( rc == INVALID_INDEX )
        {
            faceTopo->bcManager->bcRecord->bcType.push_back( this->face_solver.faceBcType[ iFace ] );
            faceTopo->bcManager->bcRecord->bcNameId.push_back( this->face_solver.faceBcKey[ iFace ] );
            ++ nBFaces;
        }
    }

    this->point_factory.InitLocalToGlobal();

}

std::unique_ptr< UnsGrid > GridElem::GenerateCalcGrid( int gridId )
{
    auto grid = ONEFLOW::CreateUnsGridUnique();
    grid->level = 0;
    grid->id = gridId;
    grid->localId = gridId;
    grid->type = UMESH;
    grid->volBcType = this->GetVolBcType();

    this->GenerateCalcGrid( *grid );
    return grid;
}

void GridElem::GenerateCalcGrid( UnsGrid & grid )
{
    grid.nCells = this->elem_feature.eTypes.size();
    grid.cellMesh->cellTopo.eTypes = this->elem_feature.eTypes;
    std::cout << "   nCells = " << grid.nCells << std::endl;

    int nNodes = this->point_factory.localToGlobal.size();
    grid.nodeMesh->CreateNodes(nNodes);
    grid.nNodes = nNodes;

    for (int iNode = 0; iNode < nNodes; ++iNode)
    {
        int globalId = this->point_factory.localToGlobal[iNode];

        Real x, y, z;
        this->point_factory.GetPoint(globalId, x, y, z);

        grid.nodeMesh->xN[iNode] = x;
        grid.nodeMesh->yN[iNode] = y;
        grid.nodeMesh->zN[iNode] = z;
    }

    this->CalcBoundaryType( grid );
    this->ReorderLink( grid );
    std::cout << "\n-->All the computing information is ready\n";
}

void GridElem::CalcBoundaryType( UnsGrid & grid )
{
    std::cout << "\n-->Set boundary condition......\n";
    grid.faceTopo = std::move( this->face_solver.faceTopo );
    grid.faceTopo->grid = &grid;
    grid.faceMesh->faceTopo = grid.faceTopo.get();
    int nFaces = grid.faceTopo->faces.size();
    std::cout << " nFaces = " << nFaces << "\n";
     
    BcRecord * bcRecord = grid.faceTopo->bcManager->bcRecord.get();
    int nBFaces = bcRecord->bcType.size();

    grid.nBFaces = nBFaces;

    std::cout << " nBFaces = " << nBFaces << "\n";

    BcTypeMap bcTypeMap;
    bcTypeMap.Init();

    IntField cgnsBcArray = bcRecord->bcType;

    IntSet originalBcSet, finalBcSet;
    int iCount = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int cgnsBcType = bcRecord->bcType[ iFace ];
        int bcNameId = bcRecord->bcNameId[ iFace ];
        int bcType = bcTypeMap.Cgns2OneFlow( cgnsBcType );

        bcRecord->bcType[ iCount ] = bcType;

        originalBcSet.insert( cgnsBcType );
        finalBcSet.insert( bcType );
        ++ iCount;
    }

    IntField nBFaceSub;

    for ( IntSet::iterator iter = originalBcSet.begin(); iter != originalBcSet.end(); ++ iter )
    {
        int iCount = 0;
        for ( int iFace = 0; iFace < nBFaces; ++ iFace )
        {
            int cgnsBcType = cgnsBcArray[ iFace ];
            if ( cgnsBcType == * iter )
            {
                ++ iCount;
            }
        }
        nBFaceSub.push_back( iCount );
    }

    std::cout << " Original Boundary Condition Number is " << originalBcSet.size() << std::endl;
    iCount = 0;
    for ( IntSet::iterator iter = originalBcSet.begin(); iter != originalBcSet.end(); ++ iter )
    {
        int oriBcType = * iter;
        std::cout << " Boundary Type = " << std::setiosflags( std::ios::right ) << std::setw( 3 ) << oriBcType;
        std::cout << "  Name = " << std::setw( 23 ) << GetCgnsBcName( oriBcType );
        std::cout << " Face = " << std::setw( 6 ) << nBFaceSub[ iCount ] << std::endl;
        ++ iCount;
    }
    std::cout << std::endl;
    std::cout << " Final Boundary Condition Number is " << finalBcSet.size() << std::endl;
    std::cout << " Boundary Type : ";
    for ( IntSet::iterator iter = finalBcSet.begin(); iter != finalBcSet.end(); ++ iter )
    {
        std::cout << * iter << " ";
    }
    std::cout << std::endl;
}

void GridElem::ReorderLink( UnsGrid & grid )
{
    FaceTopo * faceTopo = grid.faceTopo.get();

    int nFaces = faceTopo->fTypes.size();
    grid.nFaces = nFaces;

    IntField f1map( nFaces ), f2map( nFaces );
    int iCount = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rc = faceTopo->rCells[ iFace ];
        if ( rc == INVALID_INDEX )
        {
            f1map[ iFace ] = iCount;
            f2map[ iCount ] = iFace;
            ++ iCount;
        }
    }

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rc = faceTopo->rCells[ iFace ];
        if ( rc != INVALID_INDEX )
        {
            f1map[ iFace ] = iCount;
            f2map[ iCount ] = iFace;
            ++ iCount;
        }
    }
    faceTopo->facesNew.resize( nFaces );
    faceTopo->lCellsNew.resize( nFaces );
    faceTopo->rCellsNew.resize( nFaces );
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int jFace = f2map[ iFace ];
        faceTopo->facesNew[ iFace ] = faceTopo->faces[ jFace ];
        faceTopo->lCellsNew[ iFace ] = faceTopo->lCells[ jFace ];
        faceTopo->rCellsNew[ iFace ] = faceTopo->rCells[ jFace ];
    }
    faceTopo->faces = faceTopo->facesNew;
    faceTopo->lCells = faceTopo->lCellsNew;
    faceTopo->rCells = faceTopo->rCellsNew;
}

ZgridElem::ZgridElem( CgnsZbase * cgnsZbase )
    : cgnsZbase( cgnsZbase )
{
}

ZgridElem::~ZgridElem() = default;

CgnsZbase * ZgridElem::GetCgnsZbase() const
{
    return this->cgnsZbase;
}

void ZgridElem::RebindCgnsZbase( CgnsZbase * cgnsZbase ) noexcept
{
    this->cgnsZbase = cgnsZbase;
}

HXVector< std::unique_ptr< GridElem > > ZgridElem::CreateGridElements( bool multiBlock ) const
{
    HXVector< std::unique_ptr< GridElem > > data;

    if ( ! multiBlock )
    {
        HXVector< CgnsZone * > zoneViews;

        const int nOriZone = cgnsZbase->GetNZones();

        for ( int iZone = 0; iZone < nOriZone; ++ iZone )
        {
            zoneViews.push_back( cgnsZbase->GetCgnsZone( iZone ) );
        }

        const int nGridElems = 1;

        for ( int iGridElem = 0; iGridElem < nGridElems; ++ iGridElem )
        {
            data.push_back( std::make_unique< GridElem >( std::move( zoneViews ) ) );
        }
    }
    else
    {
        const int nZones = cgnsZbase->GetNZones();

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            HXVector< CgnsZone * > zoneViews;
            zoneViews.push_back( cgnsZbase->GetCgnsZone( iZone ) );

            data.push_back( std::make_unique< GridElem >( std::move( zoneViews ) ) );
        }
    }

    return data;
}

void ZgridElem::PrepareUnsCalcGrid( const HXVector< std::unique_ptr< GridElem > > & data ) const
{
    const int nGridElems = data.size();
    for ( int iGridElem = 0; iGridElem < nGridElems; ++ iGridElem )
    {
        data[ iGridElem ]->PrepareUnsCalcGrid();
    }
}

Grids ZgridElem::GenerateLocalOneFlowGrids()
{
    return this->GenerateLocalOneFlowGrids( GridConfig::FromDataBase() );
}

Grids ZgridElem::GenerateLocalOneFlowGrids( const GridConfig & config )
{
    HXVector< std::unique_ptr< GridElem > > data =
        this->CreateGridElements( config.multiBlock );
    this->PrepareUnsCalcGrid( data );

    Grids grids;
    const int nGridElems = data.size();
    grids.reserve( static_cast< std::size_t >( nGridElems ) );

    for ( int iGridElem = 0; iGridElem < nGridElems; ++ iGridElem )
    {
        grids.push_back( data[ iGridElem ]->GenerateCalcGrid( iGridElem ) );
    }

    return grids;
}


EndNameSpace
