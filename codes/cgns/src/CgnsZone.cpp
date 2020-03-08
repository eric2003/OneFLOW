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

void CgnsZone::DumpCgnsZone( Grid * grid )
{
    ONEFLOW::DumpCgnsZoneAttribute( this, grid );

    ONEFLOW::DumpCgnsGridBoundary( this, grid );

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

void CgnsZone::ReadCgnsZoneAttribute()
{
    this->ReadCgnsZoneType();

    this->ReadCgnsZoneNameAndGeneralizedDimension();

    this->SetDimension();
}

void CgnsZone::ReadCgnsZoneType()
{
    //Check the zone type
    cg_zone_type( cgnsBase->fileId, cgnsBase->baseId, this->zId, & cgnsZoneType );

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


void CgnsZone::SetDimension()
{
    this->cgnsCoor->SetDimension();
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

void CgnsZone::SetElemPosition()
{
    this->cgnsZsection->SetElemPosition();
}

void CgnsZone::ReadNumberOfCgnsSections()
{
    this->cgnsZsection->ReadNumberOfCgnsSections();
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

void CgnsZone::ReadCgnsGridBoundary()
{
    cgnsZbc->ReadCgnsGridBoundary();
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