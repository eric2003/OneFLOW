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

#include "CgnsFactory.h"
#include "CgnsGlobal.h"
#include "CgnsZbc.h"
#include "CgnsFile.h"
#include "GridTypes.h"
#include "Prj.h"
#include "Fatal.h"
#include "StringUtils.h"
#include "Su2Grid.h"
#include "GridState.h"
#include "GridMediator.h"
#include "DataBase.h"
#include "StrGrid.h"
#include "CgnsBase.h"
#include "CgnsZbase.h"
#include "CgnsZbaseUtils.h"
#include "CgnsZone.h"
#include "CgnsZoneUtils.h"
#include "CgnsSection.h"
#include "CgnsZsection.h"
#include "NodeMesh.h"
#include "PointLocator.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "HXMath.h"
#include "Dimension.h"
#include "CgnsBcBoco.h"
#include "ElementHome.h"
#include "GridDef.h"
#include "CalcGrid.h"
#include "GridElem.h"

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

// Constructor uses member initializer list and std::make_unique
CgnsFactory::CgnsFactory()
    : cgnsZbase(std::make_unique<CgnsZbase>())
{
    // Pass raw pointer to ZgridElem as it is a non-owning observer
    this->zgridElem = std::make_unique<ZgridElem>(this->cgnsZbase.get());
}

// FIX: Define destructor and move operations here.
// The compiler can now see the complete types and safely generate the 
// code to delete the unique_ptr members.
CgnsFactory::~CgnsFactory() = default;
CgnsFactory::CgnsFactory(CgnsFactory&&) noexcept = default;
CgnsFactory& CgnsFactory::operator=(CgnsFactory&&) noexcept = default;

// FIX: Exception-safe ownership transfer
void CgnsFactory::ConvertStrCgns2UnsCgnsGrid()
{
    auto unsCgnsZbase = std::make_unique<CgnsZbase>();

    // If ReadCgnsMultiBase throws, unsCgnsZbase is automatically destroyed.
    // this->cgnsZbase remains untouched and valid.
    ONEFLOW::ReadCgnsMultiBase( unsCgnsZbase.get(), this->cgnsZbase.get() );

    // Transfer ownership safely
    this->cgnsZbase = std::move(unsCgnsZbase);

    // Update the non-owning observer
    this->zgridElem->cgnsZbase = this->cgnsZbase.get();
}

void GenerateLocalOneFlowGridFromSu2Grid( Su2Grid & su2Grid, Grids & grids )
{
    // Stack allocation instead of new/delete
    CgnsFactory cgnsFactory;
    cgnsFactory.CreateSu2CgnsZone( su2Grid );

    Grids local_grids;
    cgnsFactory.zgridElem->GenerateLocalOneFlowGrid( local_grids );
    ONEFLOW::AddOneFlowGrid( grids, local_grids[ 0 ] );
}

void CgnsFactory::GenerateGrid( const std::string & caseDir )
{
    this->GenerateGrid( GridConfig::FromDataBase(), caseDir );
}

void CgnsFactory::GenerateGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    const std::string sourceCaseDir = config.sourceCaseDir.empty()
        ? caseDir
        : config.sourceCaseDir;
    this->ReadCgnsGrid( config, sourceCaseDir );

    int systemZoneType = cgnsZbase->GetSystemZoneType();
    if ( ! ( systemZoneType == CGNS_ENUMV( Unstructured ) ) )
    {
        this->ConvertStrCgns2UnsCgnsGrid();
    }

    if ( config.targetType == GridFileType::CGNS )
    {
        this->DumpUnsCgnsGrid( config, caseDir );
    }
    else
    {
        this->ProcessCgnsBases();
        this->CgnsToOneFlowGrid( config );
    }
}

void CgnsFactory::ProcessCgnsBases()
{
    this->cgnsZbase->ProcessCgnsBases();
}

void CgnsFactory::ReadCgnsGrid( const std::string & caseDir )
{
    this->ReadCgnsGrid( GridConfig::FromDataBase(), caseDir );
}

void CgnsFactory::ReadCgnsGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    // Use .get() to pass the raw pointer to legacy/global APIs
    cgns_global.cgnsbases = this->cgnsZbase.get();
    const std::string & sourceGridFile = config.sourceFile;

    std::string gridFileName;
    if ( caseDir.empty() )
    {
        gridFileName = Prj::GetPrjFileName( sourceGridFile );
    }
    else
    {
        gridFileName = Prj::GetCaseFileName( caseDir, sourceGridFile );
    }

    this->cgnsZbase->ReadCgnsGrid( gridFileName );
}

void CgnsFactory::DumpCgnsGrid( ZgridMediator & zgridMediator )
{
    cgns_global.cgnsbases = cgnsZbase.get();
    ONEFLOW::DumpCgnsGrid( cgnsZbase.get(), & zgridMediator );
}


void CgnsFactory::CommonToOneFlowGrid()
{
    this->CommonToOneFlowGrid( GridConfig::FromDataBase() );
}

void CgnsFactory::CommonToOneFlowGrid( const GridConfig & config )
{
    if ( ONEFLOW::IsUnsGrid( config.topo ) )
    {
        this->CommonToUnsGridTEST( config );
    }
    else if ( ONEFLOW::IsStrGrid( config.topo ) )
    {
        this->CommonToStrGrid();
    }
}

void CgnsFactory::CommonToStrGrid()
{
}

void CgnsFactory::DumpUnsCgnsGrid( const std::string & caseDir )
{
    this->DumpUnsCgnsGrid( GridConfig::FromDataBase(), caseDir );
}

void CgnsFactory::DumpUnsCgnsGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    const std::string & targetGridFile = config.targetFile;

    std::string targetFile;
    if ( caseDir.empty() )
    {
        targetFile = Prj::GetPrjFileName( targetGridFile );
    }
    else
    {
        targetFile = Prj::GetCaseFileName( caseDir, targetGridFile );
    }

    cgnsZbase->OpenCgnsFile( targetFile, CG_MODE_WRITE );
    cgnsZbase->DumpCgnsMultiBase();
    cgnsZbase->CloseCgnsFile();
}

void CgnsFactory::CreateCgnsZone( ZgridMediator & zgridMediator )
{
    ONEFLOW::CreateDefaultCgnsZones( cgnsZbase.get(), & zgridMediator );
}

void CgnsFactory::PrepareCgnsZone( ZgridMediator & zgridMediator )
{
    ONEFLOW::PrepareCgnsZone( cgnsZbase.get(), & zgridMediator );
}

void CgnsFactory::ReadGridAndConvertToUnsCgnsZone()
{
    this->ReadGridAndConvertToUnsCgnsZone( GridConfig::FromDataBase() );
}

void CgnsFactory::ReadGridAndConvertToUnsCgnsZone( const GridConfig & config )
{
    ZgridMediator zgridMediator;
    zgridMediator.ReadGrid( config );

    //create multi cgns zone
    this->CreateCgnsZone( zgridMediator );
    this->PrepareCgnsZone( zgridMediator );
}

void CgnsFactory::CommonToUnsGridTEST()
{
    this->CommonToUnsGridTEST( GridConfig::FromDataBase() );
}

void CgnsFactory::CommonToUnsGridTEST( const GridConfig & config )
{
    this->ReadGridAndConvertToUnsCgnsZone( config );

    this->CgnsToOneFlowGrid( config );
}

CgnsZone * CgnsFactory::CreateSu2CgnsZone( Su2Grid & su2Grid )
{
    CgnsZone * cgnsZone = this->cgnsZbase->CreateCgnsZone();

    su2Grid.FillSU2CgnsZone( *cgnsZone );

    return cgnsZone;
}

void CgnsFactory::Su2ToOneFlowGrid( Su2Grid & su2Grid )
{
    int nZones = su2Grid.nZone;
    Grids grids;

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ONEFLOW::GenerateLocalOneFlowGridFromSu2Grid( su2Grid, grids );
    }

    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void CgnsFactory::CgnsToOneFlowGrid()
{
    this->CgnsToOneFlowGrid( GridConfig::FromDataBase() );
}

void CgnsFactory::CgnsToOneFlowGrid( const GridConfig & config )
{
    if ( ! ONEFLOW::IsUnsGrid( config.topo ) ) return;

    Grids grids;

    this->zgridElem->GenerateLocalOneFlowGrid( grids );

    //The grid is processed and the grid file used for calculation is output
    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void AddOneFlowGrid( Grids & grids, Grid * grid )
{
    int iZone = grids.size() - 1;
    grids.push_back( grid );
    grid->id = iZone;
}

#endif
EndNameSpace
