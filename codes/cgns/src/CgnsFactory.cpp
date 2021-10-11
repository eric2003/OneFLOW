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
    this->zgridElem = new ZgridElem( this->cgnsZbase );
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
    cgnsZbase->OpenCgnsFile( targetFile, CG_MODE_WRITE );
    cgnsZbase->DumpCgnsMultiBase();
    cgnsZbase->CloseCgnsFile();
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

CgnsZone * CgnsFactory::CreateSu2CgnsZone( Su2Grid* su2Grid )
{
    CgnsZone * cgnsZone = this->cgnsZbase->CreateCgnsZone();

    su2Grid->FillSU2CgnsZone( cgnsZone );

    return cgnsZone;
}

void CgnsFactory::Su2ToOneFlowGrid( Su2Grid* su2Grid )
{
    int nZones = su2Grid->nZone;
    Grids grids;

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ONEFLOW::GenerateLocalOneFlowGridFromSu2Grid( su2Grid, grids );
    }

    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void CgnsFactory::CgnsToOneFlowGrid()
{
    if ( ! ONEFLOW::IsUnsGrid( grid_para.topo ) ) return;

    Grids grids;

    zgridElem->GenerateLocalOneFlowGrid( grids );

    //The grid is processed and the grid file used for calculation is output
    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void AddOneFlowGrid( Grids & grids, Grid * grid )
{
    int iZone = grids.size() - 1;
    grids.push_back( grid );
    grid->id = iZone;
}

void GenerateLocalOneFlowGridFromSu2Grid( Su2Grid* su2Grid, Grids & grids )
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->CreateSu2CgnsZone( su2Grid );

    Grids local_grids;
    cgnsFactory->zgridElem->GenerateLocalOneFlowGrid( local_grids );

    ONEFLOW::AddOneFlowGrid( grids, local_grids[ 0 ] );

    delete cgnsFactory;
}


#endif
EndNameSpace
