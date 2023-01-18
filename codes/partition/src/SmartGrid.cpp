/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

#include "SmartGrid.h"
#include "HXCgns.h"
#include "ElementHome.h"
#include "PrintDevice.h"
#include <iostream>
#include <algorithm>


BeginNameSpace( ONEFLOW )

PointAction::PointAction()
{
}

PointAction::~PointAction()
{
}

bool PointAction::FindPoint( PointAction::PointType & point, PointAction::PointMap::iterator & iter )
{
    iter = this->pointMap.find( point );
    if ( iter == this->pointMap.end() )
    {
        return false;
    }
    return true;
}

int PointAction::FindPoint( Real xm, Real ym, Real zm )
{
    PointAction::PointType pt( xm, ym, zm );
    return this->FindPointId( pt );
}

int PointAction::FindPointId( PointBasic::PointType & point )
{
    PointAction::PointMap::iterator iter = this->pointMap.find( point );
    if ( iter == this->pointMap.end() )
    {
        return -1;
    }
    return iter->second;
}

int PointAction::AddPoint( Real xm, Real ym, Real zm )
{
    PointAction::PointType pt( xm, ym, zm );
    return this->AddPoint( pt );
}

int PointAction::DeletePoint( Real xm, Real ym, Real zm )
{
    PointAction::PointType pt( xm, ym, zm );
    return this->DeletePoint( pt );
}

int PointAction::DeletePoint( PointAction::PointType & point )
{
    PointAction::PointMap::iterator iter;
    int pid = -1;
    if ( this->FindPoint( point, iter ) )
    {
        pid = iter->second;
        this->pointMap.erase( point );
        this->ModifyPointIndexAfterDelete( pid );
        int kkk = 1;
    }

    return pid;
}

void PointAction::ModifyPointIndexAfterDelete( int pid )
{
    for ( PointAction::PointMap::iterator iter = this->pointMap.begin(); iter != this->pointMap.end(); ++ iter )
    {
        if ( iter->second > pid )
        {
            iter->second -= 1;
        }
    }

    int nSize = this->pointList.size();
    int nSizeNew = nSize - 1;
    for ( int i = pid; i < nSizeNew; ++ i )
    {
        this->pointList[ i ] = this->pointList[ i + 1 ];
    }
    this->pointList.resize( nSizeNew );
}

int PointAction::AddPoint( PointAction::PointType & point )
{
    PointAction::PointMap::iterator iter;
    int pid;
    if ( this->FindPoint( point, iter ) )
    {
        pid = iter->second;
    }
    else
    {
        int index = this->pointList.size();
        point.id = index;
        pid = index;
        this->pointMap[ point ] = index;
        this->pointList.push_back( point );
    }

    return pid;
}

TopoSort::TopoSort()
{
    ;
}

TopoSort::~TopoSort()
{
    ;
}

void TopoSort::GetElementFace( UnitElement * unitElement, std::vector< int > & element, int facePos, std::vector< int > & face, int & faceType )
{
    //unit_face -> real_face mapping
    //for instance 0,1,2,3->1001,1002,1003,1004
    IntField & unit_face = unitElement->GetElementFace( facePos );
    faceType = unitElement->GetFaceType( facePos );
    int numberOfFacePoints = unit_face.size();
    face.resize( 0 );
    for ( int iFacePoint = 0; iFacePoint < numberOfFacePoints; ++ iFacePoint )
    {
        face.push_back( element[ unit_face[ iFacePoint ] ] );
    }
}

void TopoSort::AddSingleFace( UnitElement * unitElement, std::vector< int > & element, int facePos, int iCell )
{
    int faceType;
    this->GetElementFace( unitElement, element, facePos, this->real_face, faceType );

    IdTool::IDSMap::iterator iter = faceIdTool.FindIds( this->real_face, faceType );
    if ( faceIdTool.NotFind( iter ) )
    {
        int index = faceIdTool.AddData();
        this->AddNewFace( iCell, facePos, faceType );
    }
    else
    {
        int face_id = iter->second;
        this->ModifyFace( face_id, iCell, facePos );
    }
}

void TopoSort::ModifyFace( int face_id, int iCell, int face_pos )
{
    this->fBcTypes[ face_id ] = ONEFLOW::BCTypeNull; //inner bc
    this->rc[ face_id ] = iCell;
    this->rc_pos[ face_id ] = face_pos;
}

void TopoSort::AddNewFace( int iCell, int face_pos, int faceType )
{
    this->lc.push_back( iCell );
    this->rc.push_back( ONEFLOW::INVALID_INDEX );
    this->lc_pos.push_back( face_pos );
    this->rc_pos.push_back( ONEFLOW::INVALID_INDEX );
    this->fTypes.push_back( faceType );
    this->fBcTypes.push_back( ONEFLOW::INVALID_INDEX );
}

void TopoSort::AddElementFaces( std::vector< int > & element, int eType, int iCell )
{
    UnitElement * unitElement = ElementHome::GetUnitElement( eType );

    int numberOfFaceInElement = unitElement->GetElementFaceNumber();

    for ( int iLocalFace = 0; iLocalFace < numberOfFaceInElement; ++ iLocalFace )
    {
        IntField & face = unitElement->GetElementFace( iLocalFace );
        int faceType = unitElement->GetFaceType( iLocalFace );

        this->AddSingleFace( unitElement, element, iLocalFace, iCell );
    }
}

void TopoSort::ReorderFaces()
{
    std::vector<int > orderMap;
    this->CalcOrderMap( orderMap );

    this->ReOrder( this->lc, orderMap );
    this->ReOrder( this->rc, orderMap );

    this->ReOrder( this->lc_pos, orderMap );
    this->ReOrder( this->rc_pos, orderMap );

    this->ReOrder( this->fTypes, orderMap );
    this->ReOrder( this->fBcTypes, orderMap );

    this->ReOrderMapdata( this->faceIdTool, orderMap );
}

void TopoSort::CalcOrderMap( std::vector<int > &orderMap )
{
    int nFaces = this->fTypes.size();
	orderMap.resize( nFaces );

	int iBoundaryFaceCount = 0;
	int iCount = 0;
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int rc = this->rc[ iFace ];
		if ( rc == ONEFLOW::INVALID_INDEX )
		{
			orderMap[ iCount ++ ] = iFace;
			++ iBoundaryFaceCount;
		}
	}

	int nBFaces = iBoundaryFaceCount;

	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int rc = this->rc[ iFace ];
		if ( rc != ONEFLOW::INVALID_INDEX )
		{
			orderMap[ iCount ++ ] = iFace;
		}
	}
}

void TopoSort::ReOrder( std::vector< int > & varList, std::vector< int > & orderMap )
{
	std::vector< int > dataSwap = varList;
	size_t nElements = varList.size();

	for ( size_t i = 0; i < nElements; ++ i )
	{
		size_t j = orderMap[ i ];
		varList[ i ] = dataSwap[ j ];
	}
}

void TopoSort::ReOrderMapdata( IdTool & faceIdTool, std::vector< int > & orderMap )
{
    std::vector< Ids > dataSwap = faceIdTool.ids_list;
    size_t nElements = dataSwap.size();
	for ( size_t i = 0; i < nElements; ++ i )
	{
		size_t j = orderMap[ i ];
		faceIdTool.ids_list[ i ] = dataSwap[ j ];
        faceIdTool.ModifyDataIndex( dataSwap[ j ], i );
	}
}

void TopoSort::TopoPostprocess()
{
    this->ReorderFaces();
    this->ScanBcFace();
    this->SetBcGhostCell();
}

void TopoSort::ScanBcFace()
{
    int nFaces = this->fTypes.size();
    std::cout << " nFaces = " << nFaces << "\n";

	int nTraditionalBc = 0;
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int originalBcType = this->fBcTypes[ iFace ];
		if ( originalBcType == ONEFLOW::INVALID_INDEX )
		{
			++ nTraditionalBc;
		}
	}
    std::cout << " nTraditionalBc = " << nTraditionalBc << "\n";
	this->bcTypes.resize( nTraditionalBc );
}

void TopoSort::SetBcGhostCell()
{
    int nBFaces = this->bcTypes.size();

	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		this->rc[ iFace ] = iFace + nCells;
	}
}

Ids::Ids()
{
    ;
}

Ids::~Ids()
{
    ;
}

bool CompareIds::operator()( const Ids & lhs, const Ids & rhs ) const
{
    if ( lhs.type == rhs.type )
    {
        return lhs.sorted_ids < rhs.sorted_ids;
    }
    return lhs.type < rhs.type;
}

IdTool::IdTool()
{
    ;
}

IdTool::~IdTool()
{
    ;
}


IdTool::IDSMap::iterator IdTool::FindIds( const std::vector< int > & ids, int type )
{
    this->vint.type = type;
    this->vint.ids = ids;
    this->vint.sorted_ids = ids;
    std::sort( this->vint.sorted_ids.begin(), this->vint.sorted_ids.end() );
    IdTool::IDSMap::iterator iter = this->ids_map.find( this->vint );
    return iter;
}

bool IdTool::NotFind( IdTool::IDSMap::iterator & iter )
{
    return iter == this->ids_map.end();
}

int IdTool::AddData()
{
    int index = this->ids_map.size();
    this->ids_map[ this->vint ] = index;
    this->ids_list.push_back( this->vint );
    return index;
}

int IdTool::AddIds( std::vector< int > & ids, int type )
{
    IdTool::IDSMap::iterator iter = this->FindIds( ids, type );
    if ( this->NotFind( iter ) )
    {
        return this->AddData();
    }
    else
    {
        //Data already exists
        return iter->second;
    }
}

void IdTool::ModifyDataIndex( const Ids & var, int new_id )
{
    IdTool::IDSMap::iterator iter = this->FindIds( var.ids, var.type );
    if ( this->NotFind( iter ) ) return;
    iter->second = new_id;
}


TopoAction::TopoAction()
{
    this->topo_sort = new TopoSort();
}

TopoAction::~TopoAction()
{
    delete this->topo_sort;
}

void TopoAction::AddElement( int p1, int p2, int eType )
{
    std::vector< int > elem;
    elem.push_back( p1 );
    elem.push_back( p2 );

    int e_index = elementIdTool.AddIds( elem, eType );
    //std::cout << " e_index = " << e_index << "\n";
    topo_sort->AddElementFaces( elem, eType, e_index );
}

void TopoAction::CalcTopology()
{
    TopoSort topo_sort;
}


SmartGrid::SmartGrid()
{
    this->point_action = new PointAction();
    this->topo_action = new TopoAction();
}

SmartGrid::~SmartGrid()
{
    delete this->point_action;
    delete this->topo_action;
}

int SmartGrid::AddPoint( Real x, Real y, Real z )
{
    return this->point_action->AddPoint( x, y, z );
}

void SmartGrid::TestAddDeletePoints()
{
    int id1 = this->AddPoint( 0.0, 0.0, 0.0 );
    int id2 = this->AddPoint( 1.0, 0.0, 0.0 );
    int id3 = this->AddPoint( 1.0, 1.0, 1.0 );
    this->AddPoint( 2.0, 1.0, 1.0 );
    this->AddPoint( 3.0, 1.0, 1.0 );
    this->AddPoint( 4.0, 1.0, 1.0 );
    int id4 = this->point_action->DeletePoint( 1.0, 1.0, 1.0 );
    std::cout << "id1 = " << id1 << "\n";
    std::cout << "id2 = " << id2 << "\n";
    std::cout << "id3 = " << id3 << "\n";
    std::cout << "id4 = " << id4 << "\n";

    int kkk = 1;
}

void SmartGrid::Run()
{
    int ni = 41;
    int xmin = 0.0;
    int xmax = 2.0;
    this->GenerateGrid( ni, xmin, xmax );
    #ifdef ENABLE_CUDA
        InitCUDA();
    #endif
    #ifdef ENABLE_OPENMP
        #pragma omp parallel  
        std::cout << "Hello, OneFLOW OpenMP Test!\n";
    #endif
}

void SmartGrid::AddElement( int p1, int p2, int eType )
{
    this->topo_action->AddElement( p1, p2, eType );
}

void SmartGrid::GenerateGrid( int ni, Real xmin, Real xmax )
{
    Real dx = ( xmax - xmin ) / ( ni - 1 );

    std::vector< int > pid_list;
    for ( int i = 0; i < ni; ++ i )
    {
        Real xm = xmin + i * dx;
        Real ym = 0.0;
        Real zm = 0.0;

        int id = this->AddPoint( xm, ym, zm );
        pid_list.push_back( id );
    }

    int ptL = 0;
    int ptR = ni - 1;

    for ( int i = 0; i < ni - 1; ++ i )
    {
        int p1 = pid_list[ i ];
        int p2 = pid_list[ i + 1 ];

        int eType = ONEFLOW::BAR_2;

        this->AddElement( p1, p2, eType );
    }
 //   CreatBCRegion( "LeftOutFlow", ONEFLOW::BCOutflow );
	//scalarBccoL->bcName = "LeftOutFlow";
	//scalarBccoL->bcType = ONEFLOW::BCOutflow;

 //   this->AddBoundaryFace( ptl, );

    this->TopoPostprocess();

    int kkk = 1;

}

void SmartGrid::TopoPostprocess()
{
    topo_action->topo_sort->TopoPostprocess();
}

void SmartGrid::CalcTopology()
{
    this->topo_action->CalcTopology();
}

EndNameSpace
