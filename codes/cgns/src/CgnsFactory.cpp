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
#include "CgnsBcRegionProxy.h"
#include "GridPara.h"
#include "LogFile.h"
#include "Stop.h"
#include "StrUtil.h"
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
#include "CompGrid.h"
#include "GridElem.h"
#include "BgGrid.h"
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsFactory::CgnsFactory()
{
    this->cgnsMultiBase = new CgnsMultiBase();
    this->gridElemS = new GridElemS();
    this->nZone = 0;
}

CgnsFactory::~CgnsFactory()
{
    delete this->cgnsMultiBase;
    delete this->gridElemS;
}

void CgnsFactory::GenerateGrid()
{
    this->ReadCgnsGrid();
    this->CgnsToOneFlowGrid();
}

void CgnsFactory::ReadCgnsGrid()
{
    cgns_global.cgnsbases = cgnsMultiBase;
    cgnsMultiBase->ReadCgnsGrid();
}

void CgnsFactory::DumpCgnsGrid( GridMediatorS * gridMediators )
{
    cgns_global.cgnsbases = cgnsMultiBase;
    cgnsMultiBase->DumpCgnsGrid( gridMediators );
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

void CgnsFactory::CgnsToOneFlowGrid()
{
    if ( ! ONEFLOW::IsUnsGrid( grid_para.topo ) ) return;

    int systemZoneType = cgnsMultiBase->GetSystemZoneType();
    if ( ! ( systemZoneType == Unstructured ) )
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
    cgnsMultiBase->CreateDefaultCgnsZones( gridMediators );
}

void CgnsFactory::PrepareCgnsZone( GridMediatorS * gridMediators )
{
    cgnsMultiBase->PrepareCgnsZone( gridMediators );
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

void CgnsFactory::CommonToUnsGridTEST()
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

    GridMediatorS gridMediators;
    gridMediators.AddGridMediator( gridMediator );

    CgnsFactory * cgnsFactory = new CgnsFactory();
    //create multi cgns zone
    cgnsFactory->CreateCgnsZone( & gridMediators );
    cgnsFactory->PrepareCgnsZone( & gridMediators );

    //for ( int iZone = 0; iZone < nZones; ++ iZone )
    //{
    //    CgnsZone * cgnsZone = 0;

    //    Grids grid_array;

    //    if ( grid_para.multiBlock )
    //    {
    //        grid_array.push_back( gridMediator->gridVector[ iZone ] );
    //    }
    //    else
    //    {
    //        grid_array = gridMediator->gridVector;
    //    }

    //    PrepareCgnsZone( grid_array, cgnsZone );
    //    
    //    cgnsFactory->CgnsToOneFlowGrid( grids[ iZone ], iZone );
    //}

    delete cgnsFactory;

    ONEFLOW::GenerateMultiZoneCompGrids( grids );
    delete gridMediator;
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
    this->nOriZone = this->cgnsMultiBase->GetNZone();

    if ( grid_para.multiBlock == 0 )
    {
        HXVector< CgnsZone * > cgnsZones;

        for ( int iZone = 0; iZone < this->nOriZone; ++ iZone )
        {
            cgnsZones.push_back( cgnsMultiBase->GetCgnsZone( iZone ) );
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
            cgnsZones.push_back( cgnsMultiBase->GetCgnsZone( iZone ) );

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
    gridMediatorS->CreateSimple( this->nZone );
    cgnsMultiBase->CreateDefaultCgnsZones( gridMediatorS );
    delete gridMediatorS;

}

CgnsZone * CgnsFactory::CreateOneUnsCgnsZone( int cgnsZoneId )
{
    int nZone = 1;
    this->CreateDefaultZone( nZone );

    int iZone = 0;
    CgnsZone * cgnsZone = cgnsMultiBase->GetCgnsZone( iZone );
    cgnsZone->cgnsZoneType = ONEFLOW::Unstructured;
    cgnsZone->zId = cgnsZoneId;
    return cgnsZone;
}


#endif
EndNameSpace