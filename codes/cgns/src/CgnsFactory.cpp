/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "CgnsFile.h"
#include "GridPara.h"
#include "LogFile.h"
#include "Prj.h"
#include "Stop.h"
#include "StrUtil.h"
#include "Su2Grid.h"
#include "GridState.h"
#include "GridMediator.h"
#include "DataBase.h"
#include "StrGrid.h"
#include "CgnsBase.h"
#include "CgnsZbase.h"
#include "CgnsZbaseUtil.h"
#include "CgnsZone.h"
#include "CgnsZoneUtil.h"
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
#include "CalcGrid.h"
#include "GridElem.h"
#include "BgGrid.h"
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsFactory::CgnsFactory()
{
    this->cgnsZbase = new CgnsZbase();
    this->zgridElem = new ZgridElem();
    this->nZone = 0;
}

CgnsFactory::~CgnsFactory()
{
    delete this->cgnsZbase;
    delete this->zgridElem;
}

void CgnsFactory::GenerateGrid()
{
    this->ReadCgnsGrid();

    int systemZoneType = cgnsZbase->GetSystemZoneType();
    if ( ! ( systemZoneType == CGNS_ENUMV( Unstructured ) ) )
    {
        this->ConvertStrCgns2UnsCgnsGrid();
    }

    string target_filetype = grid_para.target_filetype; 
    if ( target_filetype == "cgns" )
    {
        this->DumpUnsCgnsGrid();
    }
    else
    {
        this->ProcessCgnsBases();
        this->CgnsToOneFlowGrid();
    }
}

void CgnsFactory::ProcessCgnsBases()
{
    this->cgnsZbase->ProcessCgnsBases();
}

void CgnsFactory::ReadCgnsGrid()
{
    cgns_global.cgnsbases = cgnsZbase;
    string prjFileName = ONEFLOW::GetPrjFileName( grid_para.gridFile );
    cgnsZbase->ReadCgnsGrid( prjFileName );
}

void CgnsFactory::DumpCgnsGrid( ZgridMediator * zgridMediator )
{
    cgns_global.cgnsbases = cgnsZbase;
    ONEFLOW::DumpCgnsGrid( cgnsZbase, zgridMediator );
}

void CgnsFactory::ConvertStrCgns2UnsCgnsGrid()
{
    CgnsZbase * unsCgnsZbase = new CgnsZbase();

    ONEFLOW::ReadCgnsMultiBase( unsCgnsZbase, cgnsZbase );

    delete cgnsZbase;

    cgnsZbase = unsCgnsZbase;
}

void CgnsFactory::CommonToOneFlowGrid()
{
    if ( ONEFLOW::IsUnsGrid( grid_para.topo ) )
    {
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

void CgnsFactory::DumpUnsCgnsGrid()
{
    string targetFile = ONEFLOW::GetPrjFileName( grid_para.targetFile );
    cgnsZbase->cgnsFile->OpenCgnsFile( targetFile, CG_MODE_WRITE );
    cgnsZbase->DumpCgnsMultiBase();
    cgnsZbase->cgnsFile->CloseCgnsFile();
}

void CgnsFactory::CgnsToOneFlowGrid()
{
    if ( ! ONEFLOW::IsUnsGrid( grid_para.topo ) ) return;

    this->AllocateGridElem();

    this->PrepareUnsCalcGrid();

    this->GenerateCalcGrid();

    Grids grids;
    zgridElem->GetGrids( grids );

    //The grid is processed and the grid file used for calculation is output
    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void CgnsFactory::CreateCgnsZone( ZgridMediator * zgridMediator )
{
    ONEFLOW::CreateDefaultCgnsZones( cgnsZbase, zgridMediator );
}

void CgnsFactory::PrepareCgnsZone( ZgridMediator * zgridMediator )
{
    ONEFLOW::PrepareCgnsZone( cgnsZbase, zgridMediator );
}

void CgnsFactory::ReadGridAndConvertToUnsCgnsZone()
{
    ZgridMediator zgridMediator;
    zgridMediator.ReadGrid();

    //create multi cgns zone
    this->CreateCgnsZone( & zgridMediator );
    this->PrepareCgnsZone( & zgridMediator );
}

void CgnsFactory::CommonToUnsGridTEST()
{
    this->ReadGridAndConvertToUnsCgnsZone();

    this->CgnsToOneFlowGrid();
}

void CgnsFactory::CgnsToOneFlowGrid( Grid *& grid, int zId )
{
    this->AllocateGridElem();

    this->PrepareUnsCalcGrid();

    this->GenerateCalcGrid();

    GridElem * gridElem = zgridElem->GetGridElem( 0 );
    grid = gridElem->grid;
    grid->id = zId;
}

void CgnsFactory::AllocateGridElem()
{
    this->nOriZone = this->cgnsZbase->GetNZones();

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
            zgridElem->AddGridElem( cgnsZones, iZone );
        }

    }
    else
    {
        this->nZone = this->nOriZone;

        for ( int iZone = 0; iZone < this->nZone; ++ iZone )
        {
            HXVector< CgnsZone * > cgnsZones;
            cgnsZones.push_back( cgnsZbase->GetCgnsZone( iZone ) );

            zgridElem->AddGridElem( cgnsZones, iZone );
        }
    }
}

void CgnsFactory::PrepareUnsCalcGrid()
{
    this->zgridElem->PrepareUnsCalcGrid();
}

void CgnsFactory::GenerateCalcGrid()
{
    this->zgridElem->GenerateCalcGrid();
}

void CgnsFactory::CreateDefaultZone( int nZone )
{
    ZgridMediator * zgridMediator = new ZgridMediator();
    zgridMediator->CreateSimple( nZone );
    ONEFLOW::CreateDefaultCgnsZones( cgnsZbase, zgridMediator );
    delete zgridMediator;

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

void CgnsFactory::CreateOneSu2Grid( Su2Grid* su2Grid, int iZone, Grid *& grid )
{
    int cgnsZoneId = iZone + 1;
    CgnsZone * cgnsZone = this->CreateOneUnsCgnsZone( cgnsZoneId );

    FillSU2CgnsZone( su2Grid, cgnsZone );

    this->CgnsToOneFlowGrid( grid, iZone );
}

void CgnsFactory::CreateSu2Grid( Su2Grid* su2Grid )
{
    int nZones = su2Grid->nZone;
    Grids grids( nZones );

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsFactory * cgnsFactory = new CgnsFactory();
        
        cgnsFactory->CreateOneSu2Grid(su2Grid, iZone, grids[ iZone ] );
        
        delete cgnsFactory;
    }

    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

#endif
EndNameSpace