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

#include "CgnsFactory.h"
#include "GridFactory.h"
#include "CgnsGlobal.h"
#include "CgnsZbc.h"
#include "GridPara.h"
#include "LogFile.h"
#include "Stop.h"
#include "StrUtil.h"
#include "GridState.h"
#include "GridMediator.h"
#include "DataBase.h"
#include "StrGrid.h"
#include "CgnsBase.h"
#include "CgnsZbase.h"
#include "CgnsZone.h"
#include "CgnsSection.h"
#include "CgnsZsection.h"
#include "NodeMesh.h"
#include "PointSearch.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "HXMath.h"
#include "Dimension.h"
#include "CgnsBcBoco.h"
#include "ElementHome.h"
#include "HXPointer.h"
#include "CompGrid.h"
#include "GridElem.h"
#include "BgGrid.h"
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsFactory::CgnsFactory()
{
    this->cgnsZbase = new CgnsZbase();
    this->gridElemS = new GridElemS();
    this->nZone = 0;
}

CgnsFactory::~CgnsFactory()
{
    delete this->cgnsZbase;
    delete this->gridElemS;
}

void CgnsFactory::GenerateGrid()
{
    this->ReadCgnsGrid();
    this->CgnsToOneFlowGrid();
}

void CgnsFactory::ReadCgnsGrid()
{
    cgns_global.cgnsbases = cgnsZbase;
    cgnsZbase->ReadCgnsGrid();
}

void CgnsFactory::DumpCgnsGrid( GridMediatorS * gridMediators )
{
    cgns_global.cgnsbases = cgnsZbase;
    cgnsZbase->DumpCgnsGrid( gridMediators );
}

void CgnsFactory::ConvertStrCgns2UnsCgnsGrid()
{
    CgnsZbase * unsCgnsMultiBase = new CgnsZbase();

    unsCgnsMultiBase->ConvertStrCgns2UnsCgnsGrid( cgnsZbase );

    delete cgnsZbase;

    cgnsZbase = unsCgnsMultiBase;
}

void CgnsFactory::CommonToOneFlowGrid()
{
    if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
    {
        //this->CommonToUnsGrid();
        this->CommonToUnsGridTEST();
    }
    else if ( ONEFLOW::IsStrGrid( grid_para.topo ) )
    {
        this->CommonToStrGrid();
    }
}

void CgnsFactory::CommonToStrGrid()
{
}

void CgnsFactory::CgnsToOneFlowGrid()
{
    if ( ! ONEFLOW::IsUnsGrid( grid_para.topo ) ) return;

    int systemZoneType = cgnsZbase->GetSystemZoneType();
    if ( ! ( systemZoneType == CGNS_ENUMV( Unstructured ) ) )
    {
        this->ConvertStrCgns2UnsCgnsGrid();
    }

    this->AllocateGridElem();

    this->PrepareUnsCompGrid();

    this->GenerateCompGrid();

    Grids grids;

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        GridElem * gridElem = gridElemS->GetGridElem( iZone );
        Grid * grid = gridElem->grid;
        grids.push_back( grid );
    }


    //对网格进行处理并输出计算所用的网格文件
    ONEFLOW::GenerateMultiZoneCompGrids( grids );
}

void CgnsFactory::CreateCgnsZone( GridMediatorS * gridMediators )
{
    cgnsZbase->CreateDefaultCgnsZones( gridMediators );
}

void CgnsFactory::PrepareCgnsZone( GridMediatorS * gridMediators )
{
    cgnsZbase->PrepareCgnsZone( gridMediators );
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
        int cgnsZoneId = iZone + 1;
        CgnsZone * cgnsZone = cgnsFactory->CreateOneUnsCgnsZone( cgnsZoneId );

        Grids grid_array;

        if ( grid_para.multiBlock )
        {
            grid_array.push_back( gridMediator->gridVector[ iZone ] );
        }
        else
        {
            grid_array = gridMediator->gridVector;
        }

        ONEFLOW::PrepareCgnsZone( grid_array, cgnsZone );
        
        cgnsFactory->CgnsToOneFlowGrid( grids[ iZone ], iZone );

        delete cgnsFactory;
    }

    ONEFLOW::GenerateMultiZoneCompGrids( grids );
    delete gridMediator;
}

void CgnsFactory::ReadGridAndConvertToUnsCgnsZone()
{
    GridMediatorS gridMediators;
    gridMediators.ReadGrid();

    //create multi cgns zone
    this->CreateCgnsZone( & gridMediators );
    this->PrepareCgnsZone( & gridMediators );
}

void CgnsFactory::CommonToUnsGridTEST()
{
    this->ReadGridAndConvertToUnsCgnsZone();

    this->AllocateGridElem();

    this->PrepareUnsCompGrid();

    this->GenerateCompGrid();

    Grids grids;

    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        GridElem * gridElem = gridElemS->GetGridElem( iZone );
        Grid * grid = gridElem->grid;
        grids.push_back( grid );
    }

    ONEFLOW::GenerateMultiZoneCompGrids( grids );
}

void CgnsFactory::CgnsToOneFlowGrid( Grid *& grid, int zId )
{
    this->AllocateGridElem();

    this->PrepareUnsCompGrid();

    this->GenerateCompGrid();

    GridElem * gridElem = gridElemS->GetGridElem( 0 );
    grid = gridElem->grid;
    grid->id = zId;
}

void CgnsFactory::AllocateGridElem()
{
    this->nOriZone = this->cgnsZbase->GetNZone();

    if ( grid_para.multiBlock == 0 )
    {
        HXVector< CgnsZone * > cgnsZones;

        for ( int iZone = 0; iZone < this->nOriZone; ++ iZone )
        {
            cgnsZones.push_back( cgnsZbase->GetCgnsZone( iZone ) );
        }

        this->nZone = 1;

        for ( int iZone = 0; iZone < this->nZone; ++ iZone )
        {
            gridElemS->AddGridElem( cgnsZones, iZone );
        }

    }
    else
    {
        this->nZone = this->nOriZone;

        for ( int iZone = 0; iZone < this->nZone; ++ iZone )
        {
            HXVector< CgnsZone * > cgnsZones;
            cgnsZones.push_back( cgnsZbase->GetCgnsZone( iZone ) );

            gridElemS->AddGridElem( cgnsZones, iZone );
        }
    }
}

void CgnsFactory::PrepareUnsCompGrid()
{
    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        GridElem * gridElem = gridElemS->GetGridElem( iZone );
        gridElem->PrepareUnsCompGrid();
    }
}

void CgnsFactory::GenerateCompGrid()
{
    if ( ONEFLOW::IsStrGrid( grid_para.topo ) )
    {
        this->GenerateStrCompGrid();
    }
    else if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
    {
        this->GenerateUnsCompGrid();
    }
    else
    {
        //混合网格
    }
}

void CgnsFactory::GenerateStrCompGrid()
{
}

void CgnsFactory::GenerateUnsCompGrid()
{
    for ( int iZone = 0; iZone < this->nZone; ++ iZone )
    {
        GridElem * gridElem = gridElemS->GetGridElem( iZone );
        gridElem->GenerateCompGrid();
    }
}

void CgnsFactory::CreateDefaultZone( int nZone )
{
    GridMediatorS * gridMediatorS = new GridMediatorS();
    gridMediatorS->CreateSimple( nZone );
    cgnsZbase->CreateDefaultCgnsZones( gridMediatorS );
    delete gridMediatorS;

}

CgnsZone * CgnsFactory::CreateOneUnsCgnsZone( int cgnsZoneId )
{
    int nZone = 1;
    this->CreateDefaultZone( nZone );

    int iZone = 0;
    CgnsZone * cgnsZone = cgnsZbase->GetCgnsZone( iZone );
    cgnsZone->cgnsZoneType = ONEFLOW::CGNS_ENUMV( Unstructured );
    cgnsZone->zId = cgnsZoneId;
    return cgnsZone;
}


#endif
EndNameSpace