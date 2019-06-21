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

#include "CgnsFactory.h"
#include "GridFactory.h"
#include "CgnsGlobal.h"
#include "CgnsBcRegionProxy.h"
#include "GridPara.h"
#include "LogFile.h"
#include "BasicIO.h"
#include "GridState.h"
#include "GridMediator.h"
#include "DataBase.h"
#include "StrGrid.h"
#include "CgnsBase.h"
#include "CgnsMultiBase.h"
#include "CgnsZone.h"
#include "CgnsSection.h"
#include "CgnsMultiSection.h"
#include "NodeMesh.h"
#include "PointSearch.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "HXMath.h"
#include "Dimension.h"
#include "CgnsBcRegion.h"
#include "ElementHome.h"
#include "HXPointer.h"
#include "CmpGrid.h"
#include "GridElem.h"
#include "BgGrid.h"
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsFactory::CgnsFactory()
{
    this->cgnsMultiBase = new CgnsMultiBase();
    this->nZone = 0;
}

CgnsFactory::~CgnsFactory()
{
    delete cgnsMultiBase;
}

void CgnsFactory::GenerateGrid()
{
	this->ReadCgnsGrid();
    this->ProcessGrid();
}

void CgnsFactory::ReadCgnsGrid()
{
	cgns_global.cgnsbases = cgnsMultiBase;
	cgnsMultiBase->ReadCgnsGrid();
}

void CgnsFactory::ProcessGrid()
{
	int systemZoneType = cgnsMultiBase->GetSystemZoneType();
    if ( ! ( systemZoneType == Unstructured ) &&
		 ONEFLOW::IsUnsGrid( grid_para.topo )  )
	{
		this->ConvertStrCgns2UnsCgnsGrid();
	}

	if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
	{
		this->AllocateGridElem();

        //From the standard grid information (CGNS) that we read,
        //we can get the grid information needed for the actual calculation.
        this->PrepareUnsCompGrid();
	}

    RegionNameMap::DumpRegion();

    this->AllocateCmpGrid();

    this->GenerateCmpGrid();

	//对网格进行处理并输出计算所用的网格文件
	ONEFLOW::GenerateMultiZoneCmpGrids( cmpGrids );
}

void CgnsFactory::ConvertStrCgns2UnsCgnsGrid()
{
    CgnsMultiBase * unsCgnsMultiBase = new CgnsMultiBase();

	unsCgnsMultiBase->ConvertStrCgns2UnsCgnsGrid( cgnsMultiBase );

	delete cgnsMultiBase;

	cgnsMultiBase = unsCgnsMultiBase;
}

void CgnsFactory::CommonToOneFlowGrid()
{
	if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
	{
		this->CommonToUnsGrid();
	}
	else if ( ONEFLOW::IsStrGrid( grid_para.topo ) )
	{
		this->CommonToStrGrid();
	}
}

void CgnsFactory::CommonToStrGrid()
{
}

void CgnsFactory::CommonToUnsGrid()
{
	GridMediator * gridMediator = new GridMediator();
    gridMediator->gridFile = ONEFLOW::GetDataValue< string >( "sourceGridFileName" );
    gridMediator->bcFile   = ONEFLOW::GetDataValue< string >( "sourceGridBcName" );

    gridMediator->gridType = grid_para.filetype;
	gridMediator->ReadGrid();

    int nZones = gridMediator->numberOfZones;

    if ( grid_para.multiBlock )
    {
        nZones = gridMediator->numberOfZones;
    }
    else
    {
        nZones = 1;
    }

	Grids grids( nZones );

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		CgnsFactory * cgnsFactory = new CgnsFactory();
        Grids grid_array;

        if ( grid_para.multiBlock )
        {
            grid_array.push_back( gridMediator->gridVector[ iZone ] );
        }
        else
        {
            grid_array = gridMediator->gridVector;
        }

		cgnsFactory->PrepareSection( grid_array, iZone );

		cgnsFactory->CgnsStr2Uns( grids[ iZone ], iZone );

		delete cgnsFactory;
	}

	ONEFLOW::GenerateMultiZoneCmpGrids( grids );
	delete gridMediator;
}

void CgnsFactory::CgnsStr2Uns( Grid *& grid, int zId )
{
    this->AllocateGridElem();

	this->AllocateCmpGrid();

	this->PrepareUnsCompGrid();

    RegionNameMap::DumpRegion();

    this->GenerateCmpGrid();

    grid = this->cmpGrids[ 0 ];

    grid->id = zId;

    this->DeAllocateGridElem();

}

void CgnsFactory::AllocateGridElem()
{
	this->nOriZone = this->cgnsMultiBase->nTZones;

    if ( grid_para.multiBlock == 0 )
    {
        HXVector< CgnsZone * > cgnsZones;

	    for ( int iZone = 0; iZone < this->nOriZone; ++ iZone )
	    {
			cgnsZones.push_back( cgnsMultiBase->GetZone( iZone ) );
        }

        this->nZone = 1;
		gridElems.resize( this->nZone );

	    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
	    {
            gridElems[ iZone ] = new GridElem( cgnsZones );
	    }

    }
    else
    {
		this->nZone = this->nOriZone;

		gridElems.resize( this->nZone );

	    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
	    {
            HXVector< CgnsZone * > cgnsZones;
            cgnsZones.push_back( cgnsMultiBase->GetZone( iZone ) );
            gridElems[ iZone ] = new GridElem( cgnsZones );
	    }
    }
}

void CgnsFactory::DeAllocateGridElem()
{
	for ( int iZone = 0; iZone < gridElems.size(); ++ iZone )
	{
        delete gridElems[ iZone ];
	}
}

void CgnsFactory::PrepareUnsCompGrid()
{
	for ( int iZone = 0; iZone < this->nZone; ++ iZone )
	{
        GridElem * ge = gridElems[ iZone ];
        ge->PrepareUnsCompGrid();
	}
}

void CgnsFactory::AllocateCmpGrid()
{
	if ( ONEFLOW::IsStrGrid( grid_para.topo ) )
	{
        this->nOriZone = this->cgnsMultiBase->nTZones;
		this->nZone    = this->nOriZone;
	}

    this->cmpGrids.resize( this->nZone );

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        int cgnsZoneType = this->cgnsMultiBase->GetZone( iZone )->cgnsZoneType;
        int gridType = Cgns2OneFlowZoneType( cgnsZoneType );
        Grid * grid = ONEFLOW::CreateGrid( gridType );
        grid->level = 0;
        grid->id = iZone;
        grid->localId = iZone;
        grid->type = gridType;
        grid->volBcType = this->cgnsMultiBase->volBcType;
        this->cmpGrids[ iZone ] = grid;
    }
}

void CgnsFactory::GenerateCmpGrid()
{
	if ( ONEFLOW::IsStrGrid( grid_para.topo ) )
	{
		this->GenerateStrCmpGrid();
	}
	else if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
	{
		this->GenerateUnsCmpGrid();
	}
	else
	{
		//混合网格
	}
}

void CgnsFactory::GenerateStrCmpGrid()
{
}

void CgnsFactory::GenerateUnsCmpGrid()
{
    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
		GridElem * ge = gridElems[ iZone ];

		ge->GenerateCmpGrid( cmpGrids[ iZone ] );
    }
}

void CgnsFactory::MergeToSingleZone( Grids & grids, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell )
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
        ComputeUnsId( grid, point_search, & unsId );
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

void CgnsFactory::MergeSU2ToSingleZone( Su2Grid* su2Grid, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell )
{
}

void CgnsFactory::PrepareSectionBasic( Grids & grids, CgnsZone * cgnsZone )
{
    NodeMesh * nodeMesh = cgnsZone->nodeMesh;

    int nNode, nCell;

    HXVector< Int3D * > unsIdList;

	this->MergeToSingleZone( grids, unsIdList, nodeMesh, nNode, nCell );

    cgnsZone->nNode = nNode;
    cgnsZone->nCell = nCell;
    cgnsZone->cgnsZoneType = Unstructured;

    this->FillSection( grids, unsIdList );

	cgnsZone->ConvertToInnerDataStandard();

    ONEFLOW::DeletePointer( unsIdList );
}

void CgnsFactory::PrepareSU2SectionBasic( Su2Grid* su2Grid, CgnsZone * cgnsZone )
{
}

CgnsZone * CgnsFactory::GetCreateZone( int cgnsZoneId )
{
    this->nZone = 1;

    cgnsMultiBase->Create( this->nZone );

	int iZone = 0;
	CgnsZone * cgnsZone = cgnsMultiBase->GetZone( iZone );
    cgnsZone->zId = cgnsZoneId;
    return cgnsZone;
}

void CgnsFactory::PrepareSection( Grids & grids, int cgnsZoneId )
{
    CgnsZone * cgnsZone = GetCreateZone( cgnsZoneId );
    PrepareSectionBasic( grids, cgnsZone );
}

void CgnsFactory::PrepareSU2Section( Su2Grid* su2Grid, int cgnsZoneId )
{
    CgnsZone * cgnsZone = GetCreateZone( cgnsZoneId );
    PrepareSU2SectionBasic( su2Grid, cgnsZone );
}

void CgnsFactory::FillSection( Grids & grids, HXVector< Int3D * > & unsIdList )
{
    int nTBcRegion = 0;

    int nTCell = 0;
    int nBFace = 0;

    for ( int iZone = 0; iZone < grids.size(); ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        Int3D & unsId = * unsIdList[ iZone ];

        nTCell += grid->ComputeNumberOfCell();

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        size_t nBcRegions = bcRegionGroup->regions->size();

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            BcRegion * bcRegion = ( * bcRegionGroup->regions )[ ir ];
            if ( BC::IsPoleBc( bcRegion->bcType ) ) continue;
			//if ( BC::IsNotNormalBc( bcRegion->bcType ) ) continue;
			
            nBFace += bcRegion->ComputeRegionCells();
            nTBcRegion ++;
        }
    }

	cout << " nBFace = " << nBFace << "\n";

    int iZone = 0;

	CgnsZone * cgnsZone = cgnsMultiBase->GetZone( iZone );
	
    cgnsZone->nCell = nTCell;

    cgnsZone->multiSection->nSection = 2;
	cgnsZone->multiSection->Create();

    cgnsZone->multiSection->cgnsSections[ 0 ]->startId = 1;
    cgnsZone->multiSection->cgnsSections[ 0 ]->endId   = nTCell;

    cgnsZone->multiSection->cgnsSections[ 1 ]->startId = nTCell + 1;
    cgnsZone->multiSection->cgnsSections[ 1 ]->endId   = nTCell + 1 + nBFace;

    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        cgnsZone->multiSection->cgnsSections[ 0 ]->eType = HEXA_8;
        cgnsZone->multiSection->cgnsSections[ 1 ]->eType = QUAD_4;
    }
    else
    {
        cgnsZone->multiSection->cgnsSections[ 0 ]->eType = QUAD_4;
        cgnsZone->multiSection->cgnsSections[ 1 ]->eType = BAR_2;
    }

    cgnsZone->multiSection->CreateConnList();

    CgnsBcRegionProxy * bcRegionProxy = cgnsZone->bcRegionProxy;
    bcRegionProxy->nOrdinaryBcRegion = nTBcRegion;
    bcRegionProxy->n1To1 = 0;
    bcRegionProxy->nConn = 0;
    bcRegionProxy->CreateCgnsBcRegion();

    CgnsSection * secV = cgnsZone->multiSection->cgnsSections[ 0 ];
    CgnsSection * secB = cgnsZone->multiSection->cgnsSections[ 1 ];

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
        
        IJKRange::Compute( ni, nj, nk, 0, -1 );
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

            if ( BC::IsPoleBc( bcRegion->bcType ) ) continue;
			//if ( BC::IsNotNormalBc( bcRegion->bcType ) ) continue;
            int nRegionCell = bcRegion->ComputeRegionCells();

            CgnsBcRegion * cgnsBcRegion = bcRegionProxy->cgnsBcRegions[ irc ];
            
            cgnsBcRegion->gridLocation = CellCenter;
            cgnsBcRegion->nElements    = 2;
            cgnsBcRegion->bcType       = static_cast< BCType_t >( bcTypeMap->OneFlow2Cgns( bcRegion->bcType ) );
            cgnsBcRegion->pointSetType = PointRange;
            cgnsBcRegion->CreateCgnsBcConn();
            cgnsBcRegion->connList[ 0 ] = eIdPos + 1;
            cgnsBcRegion->connList[ 1 ] = eIdPos + nRegionCell;
            string bcName = GetCgnsBcName( cgnsBcRegion->bcType );
            cgnsBcRegion->name = AddString( bcName, ir );

            eIdPos += nRegionCell;

            ONEFLOW::SetUnsBcConn( bcRegion, bConnList, pos, unsId );

            irc ++;
        }
    }

    delete bcTypeMap;
}

void CgnsFactory::FillSU2Section( Su2Grid* su2Grid )
{
}

void ComputeUnsId( StrGrid * grid, PointSearch * pointSearch, Int3D * unsId )
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
				//	int width = 8;
				//	cout << " id = " << pointIndex;
				//	cout << setw( width ) << xm;
				//	cout << setw( width ) << ym;
				//	cout << setw( width ) << zm;
				//	cout << "\n";
				//}
				

                ( * unsId )( i, j, k ) = pointIndex;
            }
        }
    }
}

int OneFlow2CgnsZoneType( int zoneType )
{
    if ( zoneType == UMESH )
    {
        return Unstructured;
    }
    else
    {
        return Structured;
    }
}

int Cgns2OneFlowZoneType( int zoneType )
{
    if ( zoneType == Unstructured )
    {
        return UMESH;
    }
    else
    {
        return SMESH;
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


#endif
EndNameSpace