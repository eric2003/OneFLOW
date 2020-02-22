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
#include "FieldImp.h"
#include "FieldAlloc.h"
#include "FieldBase.h"
#include "FieldWrap.h"
#include "UsdPara.h"
#include "Grid.h"
#include "GridState.h"
#include "DataBase.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "InterFace.h"
#include "FaceTopo.h"
#include "UNsCom.h"

BeginNameSpace( ONEFLOW )

FieldProperty::FieldProperty()
{
}

FieldProperty::~FieldProperty()
{
}

void FieldProperty::AddField( const string & fieldName, int nEqu )
{
    this->data[ fieldName ] = nEqu;
}

int FieldProperty::GetNEqu( const string & fileName )
{
    std::map< string, int >::iterator iter;
    iter = this->data.find( fileName );
    if ( iter != this->data.end() )
    {
        return iter->second;
    }
    return -1;
}

IFieldProperty::IFieldProperty()
{
}

IFieldProperty::~IFieldProperty()
{
}

void IFieldProperty::AllocateInterfaceField( int nIFace, DataStorage * dataStorage )
{
    if ( nIFace <= 0 ) return;

    for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField( dataStorage, nTEqu, nIFace, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, iter->first );
        ONEFLOW::ZeroField( field, nTEqu, nIFace );
    }
}

void IFieldProperty::UploadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
            int nEqu = iter->second;

            if ( ZoneState::zid == 0 && iter->first == "gama" )
            {
                int kkk = 1;
            }

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );
            ONEFLOW::UploadInterfaceValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void IFieldProperty::DownloadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::DownloadInterfaceValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void IFieldProperty::UploadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::UploadOversetValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void IFieldProperty::DownloadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::DownloadOversetValue( grid, targetField, iter->first, nEqu );
        }
    }
}

void IFieldProperty::DeAllocateInterfaceField( DataStorage * dataStorage )
{
    for ( std::map< string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
    {
    }
}

std::map< string, int > GFieldProperty::data;

GFieldProperty::GFieldProperty()
{
}

GFieldProperty::~GFieldProperty()
{
}

void GFieldProperty::AddField( const string & fileName, int nEqu )
{
    GFieldProperty::data[ fileName ] = nEqu;
}

int GFieldProperty::GetNEqu( const string & fileName )
{
    std::map< string, int >::iterator iter;
    iter = GFieldProperty::data.find( fileName );
    if ( iter != GFieldProperty::data.end() )
    {
        return iter->second;
    }
    return -1;
}

FieldPropertyData::FieldPropertyData()
{
    this->bcField    = new FieldProperty();
    this->faceField  = new FieldProperty();
    this->innerField = new FieldProperty();
}

FieldPropertyData::~FieldPropertyData()
{
    delete this->bcField;
    delete this->faceField;
    delete this->innerField;
}

FieldManager::FieldManager()
{
    iFieldProperty = new IFieldProperty();
    usdPara   = new UsdPara();

    strManager  = new FieldPropertyData();
    unsManager  = new FieldPropertyData();
    commManager = new FieldPropertyData();
}

FieldManager::~FieldManager()
{
    delete iFieldProperty;
    delete usdPara;

    delete strManager;
    delete unsManager;
    delete commManager;
}

void FieldManager::SetField( const string & fieldName, Real value )
{
    int nTEqu = this->commManager->innerField->GetNEqu( fieldName );

    FieldHome::SetField( fieldName, value );
}

void FieldManager::AddFaceField( const string & fieldName, int nEqu )
{
    GFieldProperty::AddField( fieldName, nEqu );
    this->commManager->faceField->AddField( fieldName, nEqu );
}

void FieldManager::AddInnerField( const string & fieldName, int nEqu )
{
    GFieldProperty::AddField( fieldName, nEqu );
    this->iFieldProperty->AddField( fieldName, nEqu );
    this->commManager->innerField->AddField( fieldName, nEqu );
}

void FieldManager::AddBcField( const string & fieldName, int nEqu )
{
    GFieldProperty::AddField( fieldName, nEqu );
    this->commManager->bcField->AddField( fieldName, nEqu );
}

void FieldManager::AddInnerField( const string & fieldName, int nEqu, int type )
{
    if ( type == 2 )
    {
        this->AddInnerField( fieldName, nEqu );
    }
    else if ( type == 0 )
    {
        unsManager->innerField->AddField( fieldName, nEqu );
    }
    else
    {
        strManager->innerField->AddField( fieldName, nEqu );
    }
}

void FieldManager::AddFaceField( const string & fieldName, int nEqu, int type )
{
    if ( type == 2 )
    {
        this->AddFaceField( fieldName, nEqu );
    }
    else if ( type == 0 )
    {
        unsManager->faceField->AddField( fieldName, nEqu );
    }
    else
    {
        strManager->faceField->AddField( fieldName, nEqu );
    }
}

void FieldManager::AddBcField( const string & fieldName, int nEqu, int type )
{
    if ( type == 2 )
    {
        this->AddBcField( fieldName, nEqu );
    }
    else if ( type == 0 )
    {
        unsManager->bcField->AddField( fieldName, nEqu );
    }
    else
    {
        strManager->bcField->AddField( fieldName, nEqu );
    }
}

void FieldManager::AllocateInnerAndBcField()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        this->AllocateInnerAndBcField( grid, this->commManager );
        this->AllocateInnerAndBcField( grid, this->unsManager );
    }
}

void FieldManager::AllocateInnerAndBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData )
{
    this->AllocateInnerField( grid, fieldPropertyData );
    this->AllocateFaceField( grid, fieldPropertyData );
    this->AllocateBcField( grid,fieldPropertyData );
}

void FieldManager::AllocateInnerField( UnsGrid * grid, FieldPropertyData * fieldPropertyData )
{
    int nTCell = grid->nCell + grid->nBFace;

    map< string, int > & data = fieldPropertyData->innerField->data;

    for ( std::map< string, int >::iterator iter = data.begin(); iter != data.end(); ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField( grid, nTEqu, nTCell, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );
        ONEFLOW::ZeroField( field, nTEqu, nTCell );
    }
}

void FieldManager::AllocateFaceField( UnsGrid * grid, FieldPropertyData * fieldPropertyData )
{
    int nFace = grid->nFace;

    map< string, int > & data = fieldPropertyData->faceField->data;

    for ( std::map< string, int >::iterator iter = data.begin(); iter != data.end(); ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField( grid, nTEqu, nFace, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

        ONEFLOW::ZeroField( field, nTEqu, nFace );
    }
}

void FieldManager::AllocateBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData )
{
    int nBFace = grid->nBFace;

    map< string, int > & data = fieldPropertyData->bcField->data;

    for ( std::map< string, int >::iterator iter = data.begin(); iter != data.end(); ++ iter )
    {
        int nTEqu = iter->second;
        ONEFLOW::CreateMRField( grid, nTEqu, nBFace, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

        ONEFLOW::ZeroField( field, nTEqu, nBFace );
    }
}

map< int, FieldManager * > * FieldFactory::data = 0;

FieldFactory::FieldFactory()
{
}

FieldFactory::~FieldFactory()
{
}

void FieldFactory::Init()
{
    if ( ! FieldFactory::data )
    {
        FieldFactory::data = new map< int, FieldManager * >();
    }
}

void FieldFactory::AddFieldManager( int sTid )
{
    map< int, FieldManager * >::iterator iter;
    FieldFactory::Init();
    iter = FieldFactory::data->find( sTid );
    if ( iter == FieldFactory::data->end() )
    {
        FieldManager * fieldManager = new FieldManager();
        ( * FieldFactory::data )[ sTid ] = fieldManager;
    }
}

FieldManager * FieldFactory::GetFieldManager( int sTid )
{
    map< int, FieldManager * >::iterator iter;
    iter = FieldFactory::data->find( sTid );
    return iter->second;
}

void FieldFactory::FreeFieldManager()
{
    map< int, FieldManager * >::iterator iter;
    for ( iter = FieldFactory::data->begin(); iter != FieldFactory::data->end(); ++ iter )
    {
        delete iter->second;
    }

    FieldFactory::data->clear();
}

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFace = interFace->nIFace;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataSend = interFace->dataSend[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        for ( int iFace = 0; iFace < nIFace; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetSId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * fieldStorage )[ iEqu ][ iFace ] = ( * field2D )[ iEqu ][ iCell ];
            }
        }
    }
}

void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFace = interFace->nIFace;
        for ( int iFace = 0; iFace < nIFace; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFace = interFace->nIFace;
        for ( int iFace = 0; iFace < nIFace; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            int iBFace = grid->interFace->i2b[ iFace ];
            int tId = grid->faceTopo->rCell[ iBFace ];

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void UploadOversetValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu )
{
}


void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu )
{
}

EndNameSpace