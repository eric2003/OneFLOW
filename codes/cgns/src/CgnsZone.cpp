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

void CgnsZone::CopyISize( CgInt * isize )
{
    for ( int i = 0; i < 9; ++ i )
    {
        this->isize[ i ] = isize[ i ];
    }
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

void CgnsZone::ReadCgnsZoneBasicInfo()
{
    this->ReadCgnsZoneType();
    this->ReadCgnsZoneNameAndGeneralizedDimension();
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

void CgnsZone::ReadCgnsGridBoundary()
{
    cgnsZbc->ReadCgnsGridBoundary();
}

void CgnsZone::ProcessPeriodicBc()
{
    ;
}

void CgnsZone::GoToZone()
{
    int fileId = this->cgnsBase->fileId;
    int baseId = this->cgnsBase->baseId;
    cg_goto( fileId, baseId,"Zone_t", this->zId, "end" );
}

void CgnsZone::GoToNode( const string & nodeName, int ith )
{
    int fileId = this->cgnsBase->fileId;
    int baseId = this->cgnsBase->baseId;
    cg_goto( fileId, baseId, "Zone_t", this->zId, nodeName.c_str(), ith, "end" );
}

void CgnsZone::GoToNode( const string & nodeNamei, int ith, const string & nodeNamej, int jth )
{
    int fileId = this->cgnsBase->fileId;
    int baseId = this->cgnsBase->baseId;
    cg_goto( fileId, baseId, "Zone_t", this->zId, nodeNamei.c_str(), ith, nodeNamej.c_str(), jth, "end" );
}

void CgnsZone::GoToNode( const string & nodeNamei, int ith, const string & nodeNamej, int jth, const string & nodeNamek, int kth )
{
    int fileId = this->cgnsBase->fileId;
    int baseId = this->cgnsBase->baseId;
    cg_goto( fileId, baseId, "Zone_t", this->zId, nodeNamei.c_str(), ith, nodeNamej.c_str(), jth, nodeNamek.c_str(), kth, "end" );
}

void CgnsZone::ReadFlowEqn()
{
    int idata[6];
    int id,ige,igm,ivm,itcm,itc,itm;
    float gamma,prandtl;
    CGNS_ENUMT(GoverningEquationsType_t) itype;
    CGNS_ENUMT(ModelType_t) mtype;

    this->GoToZone();

    cg_equationset_read( &id, &ige, &igm, &ivm, &itcm, &itc, &itm );
    cout << "Eqn dimension = " << id << "\n";

    //Read 'GoverningEquations' node
    if ( ige == 1 )
    {
        this->GoToNode( "FlowEquationSet_t", 1 );
        cg_governing_read( & itype );
        cout << " Gov eqn = " << GoverningEquationsTypeName[ itype ] << "\n";
        //Read 'DiffusionModel' node
        this->GoToNode( "FlowEquationSet_t", 1, "GoverningEquations_t", 1 );
        cg_diffusion_read( idata );
        cout << "     diffusion = ";
        cout << idata[ 0 ] << ", ";
        cout << idata[ 1 ] << ", ";
        cout << idata[ 2 ] << ", ";
        cout << idata[ 3 ] << ", ";
        cout << idata[ 4 ] << ", ";
        cout << idata[ 5 ] << "\n";
    }

    // Read gas model
    if ( igm == 1 )
    {
        this->GoToNode( "FlowEquationSet_t", 1 );
        cg_model_read( "GasModel_t", & mtype );
        cout << " Gas model type = " << ModelTypeName[ mtype ] << "\n";
        this->GoToNode( "FlowEquationSet_t", 1, "GasModel_t", 1 );
        cg_array_read_as(1,CGNS_ENUMV(RealSingle),&gamma);
        cout << "     gamma = " << gamma << "\n";
    }

    //Read turbulence closure
    if ( itc == 1 )
    {
        this->GoToNode( "FlowEquationSet_t", 1 );
        cg_model_read( "TurbulenceClosure_t", & mtype );
        cout << " Turbulence closure type = " << ModelTypeName[ mtype ] << "\n";
        this->GoToNode( "FlowEquationSet_t", 1, "TurbulenceClosure_t", 1 );
        cg_array_read_as( 1, CGNS_ENUMV(RealSingle), & prandtl );
        cout << "     turb prandtl number = " << prandtl << "\n";
    }

    if ( itm == 1 )
    {
        this->GoToNode( "FlowEquationSet_t", 1 );
        cg_model_read( "TurbulenceModel_t", & mtype );
        cout << " Turbulence model type " << ModelTypeName[ mtype ] << "\n";
    }
}

#endif
EndNameSpace