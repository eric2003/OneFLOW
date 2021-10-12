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
#include "ScalarAlloc.h"
#include "DataStorage.h"
#include "DataBase.h"
#include "FieldBase.h"
#include "ScalarZone.h"
#include "ScalarIFace.h"

BeginNameSpace( ONEFLOW )

ScalarFieldProperty::ScalarFieldProperty()
{
}

ScalarFieldProperty::~ScalarFieldProperty()
{
}

void ScalarFieldProperty::AddField( const string & fieldName, int nEqu )
{
    this->data[ fieldName ] = nEqu;
}

int ScalarFieldProperty::GetNEqu( const string & fileName )
{
    std::map< string, int >::iterator iter;
    iter = this->data.find( fileName );
    if ( iter != this->data.end() )
    {
        return iter->second;
    }
    return -1;
}

ScalarFieldAlloc::ScalarFieldAlloc()
{
}

ScalarFieldAlloc::~ScalarFieldAlloc()
{
}

ScalarFieldManager::ScalarFieldManager()
{
    this->interfaceAlloc = new ScalarFieldAlloc();
    this->inner = new ScalarFieldAlloc();
    this->faceField = new ScalarFieldAlloc();
}

ScalarFieldManager::~ScalarFieldManager()
{
    delete this->interfaceAlloc;
    delete this->inner;
    delete this->faceField;
}

void ScalarFieldManager::Init()
{
    this->interfaceAlloc->AddField( "q", 1 );
    this->inner->AddField( "q", 1 );
    this->inner->AddField( "res", 1 );

    this->faceField->AddField( "qf1", 1 );
    this->faceField->AddField( "qf2", 1 );
    this->faceField->AddField( "invflux", 1 );
}

void ScalarFieldManager::AllocateAllFields()
{
    this->AllocateInnnerField();
    this->AllocateInterfaceField();
    this->AllocateFaceField();
}

void ScalarFieldManager::AllocateInterfaceField()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    ScalarIFace * scalarIFace = grid->scalarIFace;

    int nIFaces = scalarIFace->GetNIFaces();
    //cout << " nIFaces = " << nIFaces << "\n";

    if ( nIFaces == 0 ) return;

    interfaceAlloc->AllocateField( scalarIFace->dataSend, nIFaces );
    interfaceAlloc->AllocateField( scalarIFace->dataRecv, nIFaces );
}

void ScalarFieldManager::AllocateInnnerField()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nTCells = grid->GetNTCells();

    inner->AllocateField( grid, nTCells );
}

void ScalarFieldManager::AllocateFaceField()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nFaces = grid->GetNFaces();

    faceField->AllocateField( grid, nFaces );
}

void ScalarFieldManager::UploadInterfaceField()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    ScalarIFace * scalarIFace = grid->scalarIFace;

    int nIFaces = scalarIFace->GetNIFaces();
    //cout << " nIFaces = " << nIFaces << "\n";

    if ( nIFaces == 0 ) return;

    for ( std::map< string, int >::iterator iter = this->interfaceAlloc->data.begin(); iter != this->interfaceAlloc->data.end(); ++ iter )
    {
        ONEFLOW::ScalarUploadInterfaceValue( grid, iter->first );
    }
}

void ScalarFieldManager::DownloadInterfaceField()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    ScalarIFace * scalarIFace = grid->scalarIFace;

    int nIFaces = scalarIFace->GetNIFaces();
    //cout << " nIFaces = " << nIFaces << "\n";

    if ( nIFaces == 0 ) return;

    for ( std::map< string, int >::iterator iter = this->interfaceAlloc->data.begin(); iter != this->interfaceAlloc->data.end(); ++ iter )
    {
        ONEFLOW::ScalarDownloadInterfaceValue( grid, iter->first );
    }
}

void ScalarUploadInterfaceValue( ScalarGrid * grid, const string & name )
{
    MRField * field2D = ONEFLOW::GetFieldPointer< MRField >( grid, name );

    if ( field2D == 0 ) return;

    int nEqu = field2D->GetNEqu();

    DataStorage * dataSend = grid->scalarIFace->dataSend;

    MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

    int nIFaces = grid->scalarIFace->GetNIFaces();

    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int iCell;
        grid->GetSId( iFace, iCell );

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * fieldStorage )[ iEqu ][ iFace ] = ( * field2D )[ iEqu ][ iCell ];
        }
    }
}

void ScalarDownloadInterfaceValue( ScalarGrid * grid, const string & name )
{
    MRField * field2D = ONEFLOW::GetFieldPointer< MRField >( grid, name );

    if ( field2D == 0 ) return;

    int nEqu = field2D->GetNEqu();

    DataStorage * dataRecv = grid->scalarIFace->dataRecv;

    MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

    int nIFaces = grid->scalarIFace->GetNIFaces();

    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int iCell;
        grid->GetTId( iFace, iCell );

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
        }
    }
}


EndNameSpace
