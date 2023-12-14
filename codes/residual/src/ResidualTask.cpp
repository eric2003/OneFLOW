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
#include "ResidualTask.h"
#include "ActionState.h"
#include "Zone.h"
#include "ZoneState.h"
#include "PIO.h"
#include "Stop.h"
#include "Parallel.h"
#include "SolverState.h"
#include "SolverInfo.h"
#include "DataBase.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "CellMesh.h"
#include "Iteration.h"
#include "FileIO.h"
#include "StrUtil.h"
#include "UNsCom.h"
#include <sstream>
#include <iomanip>


BeginNameSpace( ONEFLOW )

ResidualTask::ResidualTask()
{
    ;
}

ResidualTask::~ResidualTask()
{
    ;
}

void ResidualTask::Run()
{
    ActionState::dataBook = this->dataBook;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( SolverState::tid );
    data.Init( solverInfo->nEqu );

    dataList.resize( ZoneState::nLocal );
    int iCount = 0;
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        if (  ! ZoneState::IsValidZone( zId ) ) continue;
        ZoneState::zid = zId;
        dataList[ iCount ].Init( solverInfo->nEqu );
        this->CalcRes( SolverState::tid, dataList[ iCount ] );
        ++ iCount;
    }

    PostDumpResiduals();
}

void ResidualTask::CalcRes( int sTid, ResData & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );

    int nEqu = solverInfo->nEqu;

    MRField * res = ONEFLOW::GetFieldPointer< MRField >( grid, solverInfo->residualName );

    data.resave.Zero();
    data.resave.nCells =  grid->nCells;

    data.resmax.index = 0;
    data.resmax.resmax = 0;

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        for ( int cId = 0; cId < grid->nCells; ++ cId )
        {
            Real ress = ( * res )[ iEqu ][ cId ];
            if ( NotANumber( ress ) )
            {
                std::cout << " iEqu = " << iEqu << " cId = " << cId << " grid->nCells = " << grid->nCells << "\n";
                std::cout << " ress = " << ress << "\n";
                Stop( "NotANumber" );
            }
            data.resave.res[ iEqu ] += SQR( ress );
            if ( data.resmax.resmax[ iEqu ] < ABS( ress ) )
            {
                data.resmax.resmax[ iEqu ] = ABS( ress );
                data.resmax.index [ iEqu ] = cId;
            }
        }
    }

    RealField & xcc = grid->cellMesh->xcc;
    RealField & ycc = grid->cellMesh->ycc;
    RealField & zcc = grid->cellMesh->zcc;
    RealField & vol = grid->cellMesh->vol;

    data.resmax.zid = grid->id;

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int id = data.resmax.index[ iEqu ];
        data.resmax.xcc[ iEqu ] = xcc[ id ];
        data.resmax.ycc[ iEqu ] = ycc[ id ];
        data.resmax.zcc[ iEqu ] = zcc[ id ];
        data.resmax.vol[ iEqu ] = vol[ id ];
    }

}

void ResidualTask::PostDumpResiduals()
{
    size_t nEqu = this->data.resave.res.size();
    this->data.resave.CalcAver( dataList );
    this->data.resmax.CalcMax( dataList );

    if ( Parallel::pid != Parallel::serverid ) return;

    this->DumpScreen();
    this->DumpFile();
}

void ResidualTask::DumpFile()
{
    std::ostringstream oss;

    std::fstream file;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( SolverState::tid );
    std::string & fileName = solverInfo->resFileName;
    PIO::OpenPrjFile( file, fileName, std::ios_base::out | std::ios_base::app );

    if ( IsEmpty( file ) )
    {
        StringField title;
        title.push_back( "Title=\"THE RESIDUAL OF ONEFLOW\"" );
        title.push_back( "Variables=" );
        title.push_back( "\"iter\"" );
        title.push_back( "\"sub-iter\"" );
        size_t nVar = this->data.resave.res.size();
        for ( int iVar = 0; iVar < nVar; ++ iVar )
        {
            title.push_back( AddString( "\"res",  iVar + 1, "\"" ) );
        }

        for ( HXSize_t iTitle = 0; iTitle < title.size(); ++ iTitle )
        {
            oss << title[ iTitle ] << std::endl;
        }
    }

    oss << std::setiosflags( std::ios::left );
    oss << std::setprecision( 5 );
    oss << std::setiosflags( std::ios::scientific );
    oss << std::setiosflags( std::ios::showpoint );

    oss << Iteration::outerSteps << " ";
    oss << Iteration::innerSteps << " ";

    size_t nVar = this->data.resave.res.size();
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        oss << std::setw( 13 ) << this->data.resave.res[ iVar ] << " ";
    }

    oss << std::endl;

    file << oss.str();

    PIO::CloseFile( file );

}

void ResidualTask::DumpScreen()
{
    int maxId = this->data.resmax.CalcMaxId();

    std::ostringstream oss;
    if ( ( Iteration::outerSteps - 1 ) % 100 == 0 )
    {
        oss << std::endl;
        oss << "iter initer ave  max zone cell vol  nv \n";
    }

    oss << std::setiosflags( std::ios::left );
    oss << std::setprecision( 5 );
    oss << std::setiosflags( std::ios::scientific );
    oss << std::setiosflags( std::ios::showpoint );

    oss << std::setw( 7  ) << Iteration::outerSteps;
    oss << std::setw( 7  ) << Iteration::innerSteps;
    oss << std::setw( 13 ) << this->data.resave.res[ maxId ];
    oss << std::setw( 13 ) << this->data.resmax.resmax[ maxId ];
    oss << std::setw( 4  ) << this->data.resmax.zid[ maxId ] + 1;
    oss << " ";
    oss << std::setw( 6  ) << this->data.resmax.index[ maxId ];
    oss << " ";
    oss << std::resetiosflags( std::ios::scientific );
    oss << std::setprecision( 3 );
    oss << std::setw( 11 ) << this->data.resmax.vol[ maxId ];
    oss << " ";
    oss << std::setw( 3 )  << maxId + 1;
    oss << std::endl;

    std::cout << oss.str();
}


EndNameSpace
