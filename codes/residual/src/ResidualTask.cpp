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
#include "ResidualTask.h"
#include "ActionState.h"
#include "Zone.h"
#include "Parallel.h"
#include "SolverState.h"
#include "SolverInfo.h"
#include "DataBase.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "CellMesh.h"
#include "Iteration.h"
#include "AsciiFileIO.h"
#include "UNsCom.h"
#include <sstream>
#include <iomanip>
using namespace std;

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
    cout << "ResidualTask::Run\n";
    ActionState::dataBook = this->dataBook;
    cout << "ResidualTask::Run  1 \n";
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( SolverState::tid );
    data.Init( solverInfo->nEqu );
    cout << "ResidualTask::Run  2 \n";

    dataList.resize( ZoneState::nLocal );
    int iCount = 0;
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        if (  ! ZoneState::IsValidZone( zId ) ) continue;
        ZoneState::zid = zId;
        dataList[ iCount ].Init( solverInfo->nEqu );
        this->CmpRes( SolverState::tid, dataList[ iCount ] );
        ++ iCount;
    }

    cout << "ResidualTask::Run  3 \n";

    PostDumpResiduals();
    cout << "ResidualTask::Run  4 \n";
}

void ResidualTask::CmpRes( int sTid, ResData & data )
{
    cout << "ResidualTask::CmpRes 1 \n";
    UnsGrid * grid = Zone::GetUnsGrid();
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );

    cout << "ResidualTask::CmpRes 2 \n";

    int nEqu = solverInfo->nEqu;

    MRField * res = ONEFLOW::GetFieldPointer< MRField >( grid, solverInfo->residualName );

    data.resave.Zero();
    data.resave.nCell =  grid->nCell;

    data.resmax.index = 0;
    data.resmax.resmax = 0;

    cout << "ResidualTask::CmpRes 3 \n";

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        for ( int cId = 0; cId < grid->nCell; ++ cId )
        {
            Real ress = ( * res )[ iEqu ][ cId ];
            data.resave.res[ iEqu ] += SQR( ress );
            if ( data.resmax.resmax[ iEqu ] < ABS( ress ) )
            {
                data.resmax.resmax[ iEqu ] = ABS( ress );
                data.resmax.index [ iEqu ] = cId;
            }
        }
    }

    cout << "ResidualTask::CmpRes 4 \n";

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

    cout << "ResidualTask::CmpRes 5 \n";
}

void ResidualTask::PostDumpResiduals()
{
    cout << "PostDumpResiduals 1 \n";
    size_t nEqu = this->data.resave.res.size();
    cout << "PostDumpResiduals 2 \n";
    this->data.resave.CmpAver( dataList );
    this->data.resmax.CmpMax( dataList );
    cout << "PostDumpResiduals 3 \n";

    if ( Parallel::pid != Parallel::serverid ) return;
    cout << "PostDumpResiduals 4 \n";

    this->DumpScreen();
    cout << "PostDumpResiduals 5 \n";
    this->DumpFile();
    cout << "PostDumpResiduals 6 \n";
}

void ResidualTask::DumpFile()
{
    cout << "DumpFile 1 \n";
    ostringstream oss;

    fstream file;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( SolverState::tid );
    string & fileName = solverInfo->resFileName;
    cout << "DumpFile 2 \n";
    PIO::ParallelOpenPrj( file, fileName, ios_base::out | ios_base::app );
    cout << "DumpFile 3 \n";

    if ( IsEmpty( file ) )
    {
        cout << "DumpFile 4 \n";
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
        cout << "DumpFile 5 \n";

        for ( UInt iTitle = 0; iTitle < title.size(); ++ iTitle )
        {
            oss << title[ iTitle ] << endl;
        }
        cout << "DumpFile 6 \n";
    }
    cout << "DumpFile 7 \n";

    oss << setiosflags( ios::left );
    oss << setprecision( 5 );
    oss << setiosflags( ios::scientific );
    oss << setiosflags( ios::showpoint );

    oss << Iteration::outerSteps << " ";
    oss << Iteration::innerSteps << " ";

    cout << "DumpFile 8 \n";

    size_t nVar = this->data.resave.res.size();
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        oss << setw( 13 ) << this->data.resave.res[ iVar ] << " ";
    }

    cout << "DumpFile 9 \n";

    oss << endl;

    file << oss.str();

    cout << "DumpFile 10 \n";

    PIO::ParallelClose( file );

    cout << "DumpFile 11 \n";
}

void ResidualTask::DumpScreen()
{
    cout << "DumpScreen 1 \n";
    int maxId = this->data.resmax.CmpMaxId();

    cout << "DumpScreen 2 \n";

    ostringstream oss;
    if ( ( Iteration::outerSteps - 1 ) % 100 == 0 )
    {
        oss << endl;
        oss << "iter initer ave  max zone cell vol  nv \n";
    }
    cout << "DumpScreen 3 \n";
    oss << setiosflags( ios::left );
    oss << setprecision( 5 );
    oss << setiosflags( ios::scientific );
    oss << setiosflags( ios::showpoint );

    oss << setw( 7  ) << Iteration::outerSteps;
    oss << setw( 7  ) << Iteration::innerSteps;
    oss << setw( 13 ) << this->data.resave.res[ maxId ];
    oss << setw( 13 ) << this->data.resmax.resmax[ maxId ];
    oss << setw( 4  ) << this->data.resmax.zid[ maxId ] + 1;
    oss << " ";
    oss << setw( 6  ) << this->data.resmax.index[ maxId ];
    oss << " ";
    oss << resetiosflags( ios::scientific );
    oss << setprecision( 3 );
    oss << setw( 11 ) << this->data.resmax.vol[ maxId ];
    oss << " ";
    oss << setw( 3 )  << maxId + 1;
    oss << endl;
    cout << "DumpScreen 4 \n";

    cout << oss.str();
    cout << "DumpScreen 5 \n";
}


EndNameSpace