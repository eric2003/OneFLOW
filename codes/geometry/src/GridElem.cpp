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

GridElem::GridElem( HXVector< std::reference_wrapper< CgnsZone > > zoneViews )
    : zoneViews( std::move( zoneViews ) ),
      minLen( LARGE ),
      maxLen( -LARGE )
{
}

GridElem::~GridElem() = default;

CgnsZone & GridElem::GetCgnsZone( int iZone )
{
    return this->zoneViews[ iZone ].get();
}

const CgnsZone & GridElem::GetCgnsZone( int iZone ) const
{
    return this->zoneViews[ iZone ].get();
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

    const bool hasPolygon = this->GetCgnsZone( 0 ).RequireCgnsZsection().HasPolygonSection();

    for ( int iZone = 1; iZone < this->GetNZones(); ++ iZone )
    {
        const bool zoneHasPolygon =
            this->GetCgnsZone( iZone ).RequireCgnsZsection().HasPolygonSection();

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

    const int volBcType = this->GetCgnsZone( 0 ).GetVolBcType();

    for ( int iZone = 1; iZone < this->GetNZones(); ++ iZone )
    {
        if ( this->GetCgnsZone( iZone ).GetVolBcType() != volBcType )
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
this->InitCgnsElements();
this->elemFeature.ScanElements( this->faceSolver );
this->ScanBcFace();
    this->GenerateCalcElement();
}

void GridElem::PrepareUnsCalcGridPolyhedron()
{
this->ScanPolygonFace();
    this->ScanBcFace();
    this->GenerateCalcElement();
}

void GridElem::ScanPolygonFace()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone & cgnsZone = this->GetCgnsZone( iZone );

        cgnsZone.ConstructCgnsGridPoints( &this->pointFactory );

        //Scan NGON_n PolygonFace
        const int nSections = cgnsZone.RequireCgnsZsection().GetNSections();
        for ( int iSection = 0; iSection < nSections; ++ iSection )
        {
            CgnsSection & cgnsSection = cgnsZone.RequireCgnsZsection().GetCgnsSection( iSection );
            if ( cgnsSection.eType != NGON_n ) continue;
            this->faceSolver.ScanPolygonFace( cgnsSection );
        }
        //Scan NFACE_n PolyhedronElement
        for ( int iSection = 0; iSection < nSections; ++ iSection )
        {
            CgnsSection & cgnsSection = cgnsZone.RequireCgnsZsection().GetCgnsSection( iSection );
            if ( cgnsSection.eType != NFACE_n ) continue;
            this->faceSolver.ScanPolyhedronElement( cgnsSection );
            this->SetPolyhedronElementType( cgnsSection );
        }

        int nFaces = this->faceSolver.GetFaceTopo().GetFaces().size();

    }
}

void GridElem::SetPolyhedronElementType( CgnsSection & cgnsSection )
{
    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        int e_type = cgnsSection.eTypeList[ iElem ];

        this->elemFeature.eTypes.push_back( e_type );
    }
}

void GridElem::InitCgnsElements()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone & cgnsZone = this->GetCgnsZone( iZone );
        
        cgnsZone.ConstructCgnsGridPoints( &this->pointFactory );
        cgnsZone.SetElementTypeAndNode( this->elemFeature );
    }
}

void GridElem::ScanBcFace()
{
    int nZone = this->GetNZones();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        CgnsZone & cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone.ScanBcFace( this->faceSolver );
    }

    this->faceSolver.ScanInterfaceBc();
}

void GridElem::GenerateCalcElement()
{
    FaceTopo & faceTopo = this->faceSolver.GetFaceTopo();
    this->faceSolver.AppendBoundaryBcToRecord( faceTopo.GetBcRecord() );

    this->pointFactory.InitLocalToGlobal();
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
    grid.nCells = this->elemFeature.eTypes.size();
    grid.GetCellMesh().GetCellTopo().eTypes = this->elemFeature.eTypes;

    int nNodes = this->pointFactory.localToGlobal.size();
    grid.nodeMesh->CreateNodes(nNodes);
    grid.nNodes = nNodes;

    for (int iNode = 0; iNode < nNodes; ++iNode)
    {
        int globalId = this->pointFactory.localToGlobal[iNode];

        Real x, y, z;
        this->pointFactory.GetPoint(globalId, x, y, z);

        grid.nodeMesh->SetPoint( iNode, x, y, z );
    }

    this->CalcBoundaryType( grid );
    this->ReorderLink( grid );
}

void GridElem::CalcBoundaryType( UnsGrid & grid )
{
    grid.SetFaceTopo( this->faceSolver.ReleaseFaceTopo() );

    BcRecord & bcRecord = grid.GetFaceTopo().GetBcRecord();
    grid.nBFaces = static_cast< int >( bcRecord.bcType.size() );

    // Map staged CGNS BC type ids to OneFLOW ids. Logging/summary used to live
    // here and mixed observation with mutation; mapping is now a pure step.
    BcTypeMap bcTypeMap;
    bcTypeMap.Init();
    bcTypeMap.MapFromCgnsInPlace( bcRecord.bcType );
}

void GridElem::ReorderLink( UnsGrid & grid )
{
    FaceTopo & faceTopo = grid.GetFaceTopo();

    int nFaces = faceTopo.GetFaceTypes().size();
    grid.nFaces = nFaces;

    faceTopo.ReorderLink();
    // BcRecord was filled in original face order of boundary faces; ReorderLink
    // preserves that relative order in the boundary prefix. Verify the contract.
    ValidateBoundaryBcLayout( faceTopo );
}

GridElemGroup::GridElemGroup( CgnsZbase & cgnsZbase )
    : cgnsZbase( cgnsZbase )
{
}

GridElemGroup::~GridElemGroup() = default;

CgnsZbase & GridElemGroup::GetCgnsZbase() const
{
    return this->cgnsZbase;
}

HXVector< std::unique_ptr< GridElem > > GridElemGroup::CreateGridElements(
    GridAssemblyMode assemblyMode ) const
{
    HXVector< std::unique_ptr< GridElem > > data;

    if ( assemblyMode == GridAssemblyMode::AggregateZones )
    {
        HXVector< std::reference_wrapper< CgnsZone > > zoneViews;

        const int nOriZone = cgnsZbase.GetNZones();

        for ( int iZone = 0; iZone < nOriZone; ++ iZone )
        {
            zoneViews.emplace_back( *cgnsZbase.GetCgnsZone( iZone ) );
        }

        const int nGridElems = 1;

        for ( int iGridElem = 0; iGridElem < nGridElems; ++ iGridElem )
        {
            data.push_back( std::make_unique< GridElem >( std::move( zoneViews ) ) );
        }
    }
    else
    {
        const int nZones = cgnsZbase.GetNZones();

        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            HXVector< std::reference_wrapper< CgnsZone > > zoneViews;
            zoneViews.emplace_back( *cgnsZbase.GetCgnsZone( iZone ) );

            data.push_back( std::make_unique< GridElem >( std::move( zoneViews ) ) );
        }
    }

    return data;
}

void GridElemGroup::PrepareUnsCalcGrid( const HXVector< std::unique_ptr< GridElem > > & data ) const
{
    const int nGridElems = data.size();
    for ( int iGridElem = 0; iGridElem < nGridElems; ++ iGridElem )
    {
        data[ iGridElem ]->PrepareUnsCalcGrid();
    }
}

Grids GridElemGroup::GenerateLocalOneFlowGrids()
{
    return this->GenerateLocalOneFlowGrids( GridConfig::FromDataBase() );
}

Grids GridElemGroup::GenerateLocalOneFlowGrids( const GridConfig & config )
{
    HXVector< std::unique_ptr< GridElem > > data =
        this->CreateGridElements( config.assemblyMode );
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
