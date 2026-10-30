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

#include "ParaFile.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Parallel.h"
#include "LogFile.h"
#include "PIO.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

BeginNameSpace( ONEFLOW )

void DumpDataBase()
{
    DataBase & dataBase = ONEFLOW::RequireGlobalDataBase();
    std::fstream file;
    std::string fileName = "/log/database.log";
    PIO::OpenPrjFile( file, fileName, std::ios_base::out );
    dataBase.RequireDataPara().DumpData( file );
    PIO::CloseFile( file );
}

void DumpDataBase( const std::string & caseDir )
{
    DataBase & dataBase = ONEFLOW::RequireGlobalDataBase();
    std::fstream file;
    PIO::OpenCaseFile( file, caseDir, "log/database.log", std::ios_base::out );
    dataBase.RequireDataPara().DumpData( file );
    PIO::CloseFile( file );
}

void BroadcastControlParameterToAllProcessors()
{
    ONEFLOW::logFile << "Broadcast Control Parameter To All Processors\n";

    ONEFLOW::HXBcast( ONEFLOW::CompressData, ONEFLOW::DecompressData, Parallel::GetServerid() );
}

void CompressData( DataBook * dataBook )
{
    DataBase & globalDataBase = ONEFLOW::RequireGlobalDataBase();

    ONEFLOW::CompressData( &globalDataBase, dataBook );
}

void DecompressData( DataBook * dataBook )
{
    DataBase & globalDataBase = ONEFLOW::RequireGlobalDataBase();
    ONEFLOW::DecompressData( &globalDataBase, dataBook );
}

void CompressData( DataBase * dataBase, DataBook * dataBook )
{
    if ( dataBase == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    const DataPara::DataMap & dataMap = dataBase->RequireDataPara().GetDataMap();

    int ndata = static_cast<int>( dataMap.size() );
    ONEFLOW::HXWrite( dataBook, ndata );

    for ( const auto & pair : dataMap )
    {
        const DataEntry * dataEntry = pair.second.get();
        ONEFLOW::HXWriteDataEntry( dataBook, dataEntry );
    }
}

void DecompressData( DataBase * dataBase, DataBook * dataBook )
{
    if ( dataBase == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    DataPara & dataPara = dataBase->RequireDataPara();
    dataBook->MoveToBegin();

    int ndata = 0;
    ONEFLOW::HXRead( dataBook, ndata );

    for ( int i = 0; i < ndata; ++ i )
    {
        auto dataEntry = ONEFLOW::HXReadDataEntry( dataBook );
        dataPara.SetDataEntry( std::move( dataEntry ) );
    }
}

EndNameSpace
