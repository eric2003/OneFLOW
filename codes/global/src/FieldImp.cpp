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
#include "FieldImp.h"
#include "FieldBase.h"
#include "FieldWrap.h"
#include "Fatal.h"
#include "UsdPara.h"
#include "Grid.h"
#include "GridState.h"
#include "DataBase.h"
#include "DataStorage.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "InterFace.h"
#include "FaceTopo.h"
#include "UNsCom.h"

BeginNameSpace( ONEFLOW )

namespace
{
    void DumpFieldProperty(
        std::ostream & output,
        const char * name,
        const FieldProperty & fieldProperty )
    {
        output
            << "  "
            << name
            << ":\n";

        if ( fieldProperty.GetData().empty() )
        {
            output << "    <empty>\n";
            return;
        }

        fieldProperty.Dump( output );
    }
}

void FieldProperty::AddField(
    const std::string & fieldName,
    int nEqu )
{
    FieldProperty::Data::iterator iter =
        this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        this->data[ fieldName ] = nEqu;
        return;
    }

    if ( iter->second != nEqu )
    {
        Fatal(
            "Conflicting field definition: "
            + fieldName );
    }
}

const FieldProperty::Data & FieldProperty::GetData() const
{
    return this->data;
}

void FieldProperty::Dump(
    std::ostream & output ) const
{
    for ( FieldProperty::Data::const_iterator iter =
        this->data.begin();
        iter != this->data.end();
        ++ iter )
    {
        output
            << "    "
            << iter->first
            << "  nEqu="
            << iter->second
            << '\n';
    }
}

void IFieldProperty::AllocateInterfaceField( int nIFaces, DataStorage * dataStorage )
{
    if ( nIFaces <= 0 ) return;

    const FieldProperty::Data & data = this->GetData();
    for ( FieldProperty::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField( dataStorage, nTEqu, nIFaces, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, iter->first );
        ONEFLOW::ZeroField( field, nTEqu, nIFaces );
    }
}

void IFieldProperty::UploadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        const FieldProperty::Data & data = this->GetData();
        for ( FieldProperty::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

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

        const FieldProperty::Data & data = this->GetData();
        for ( FieldProperty::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
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

        const FieldProperty::Data & data = this->GetData();
        for ( FieldProperty::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
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

        const FieldProperty::Data & data = this->GetData();
        for ( FieldProperty::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::DownloadOversetValue( grid, targetField, iter->first, nEqu );
        }
    }
}

FieldProperty & FieldPropertyData::GetFieldProperty(
    FieldLocation location )
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

const FieldProperty & FieldPropertyData::GetFieldProperty(
    FieldLocation location ) const
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

FieldManager::FieldManager()
    : fieldDefinitionsReady( false )
    , interfaceDefinitionsReady( false )
{
    usdPara =
        std::make_unique< UsdPara >();
}

FieldManager::~FieldManager() = default;

bool FieldManager::HasFieldDefinitions() const
{
    return this->fieldDefinitionsReady;
}

void FieldManager::MarkFieldDefinitionsReady()
{
    this->fieldDefinitionsReady = true;
}

bool FieldManager::HasInterfaceDefinitions() const
{
    return this->interfaceDefinitionsReady;
}

void FieldManager::MarkInterfaceDefinitionsReady()
{
    this->interfaceDefinitionsReady = true;
}

IFieldProperty & FieldManager::GetInterfaceFieldProperty()
{
    return this->iFieldProperty;
}

const IFieldProperty & FieldManager::GetInterfaceFieldProperty() const
{
    return this->iFieldProperty;
}

FieldPropertyData & FieldManager::GetFieldPropertyData(
    FieldCategory category )
{
    switch ( category )
    {
    case FieldCategory::Common:
        return commManager;

    case FieldCategory::Structured:
        return strManager;

    case FieldCategory::Unstructured:
        return unsManager;
    }

    Fatal( "Invalid field category" );
    return commManager;
}

const FieldPropertyData & FieldManager::GetFieldPropertyData(
    FieldCategory category ) const
{
    switch ( category )
    {
    case FieldCategory::Common:
        return commManager;

    case FieldCategory::Structured:
        return strManager;

    case FieldCategory::Unstructured:
        return unsManager;
    }

    Fatal( "Invalid field category" );
    return commManager;
}

UsdPara & FieldManager::GetUsdPara()
{
    return *this->usdPara;
}

const UsdPara & FieldManager::GetUsdPara() const
{
    return *this->usdPara;
}

void FieldManager::DumpFieldEnvironment(
    std::ostream & output ) const
{
    output
        << "========== Field Environment ==========\n\n";

    output
        << "[Common]\n";

    DumpFieldProperty(
        output,
        "Inner",
        this->commManager.GetFieldProperty(
            FieldLocation::Inner ) );

    DumpFieldProperty(
        output,
        "Face",
        this->commManager.GetFieldProperty(
            FieldLocation::Face ) );

    DumpFieldProperty(
        output,
        "Boundary",
        this->commManager.GetFieldProperty(
            FieldLocation::Boundary ) );

    output
        << "\n[Structured]\n";

    DumpFieldProperty(
        output,
        "Inner",
        this->strManager.GetFieldProperty(
            FieldLocation::Inner ) );

    DumpFieldProperty(
        output,
        "Face",
        this->strManager.GetFieldProperty(
            FieldLocation::Face ) );

    DumpFieldProperty(
        output,
        "Boundary",
        this->strManager.GetFieldProperty(
            FieldLocation::Boundary ) );

    output
        << "\n[Unstructured]\n";

    DumpFieldProperty(
        output,
        "Inner",
        this->unsManager.GetFieldProperty(
            FieldLocation::Inner ) );

    DumpFieldProperty(
        output,
        "Face",
        this->unsManager.GetFieldProperty(
            FieldLocation::Face ) );

    DumpFieldProperty(
        output,
        "Boundary",
        this->unsManager.GetFieldProperty(
            FieldLocation::Boundary ) );

    output
        << "\n[Interface Storage]\n";

    if ( this->iFieldProperty.GetData().empty() )
    {
        output << "    <empty>\n";
    }
    else
    {
        this->iFieldProperty.Dump( output );
    }

    output
        << "\n========================================\n";
}

void FieldManager::SetField( const std::string & fieldName, Real value )
{
    FieldHome::SetField( fieldName, value );
}

void FieldManager::AddField(
    const std::string & fieldName,
    int nEqu,
    FieldCategory category,
    FieldLocation location )
{
    FieldProperty & fieldProperty =
        this->GetFieldPropertyData(
            category ).GetFieldProperty(
                location );

    fieldProperty.AddField(
        fieldName,
        nEqu );
}

void FieldManager::AddInterfaceField(
    const std::string & fieldName )
{
    const FieldProperty::Data & data =
        this->commManager.GetFieldProperty(
            FieldLocation::Inner ).GetData();

    FieldProperty::Data::const_iterator iter =
        data.find( fieldName );

    if ( iter == data.end() )
    {
        Fatal(
            "Interface field is not defined as a common inner field: "
            + fieldName );
    }

    this->iFieldProperty.AddField(
        fieldName,
        iter->second );
}

std::map< int, std::unique_ptr< FieldManager > > FieldFactory::data;

void FieldFactory::AddFieldManager( int solverType )
{
    std::map< int, std::unique_ptr< FieldManager > >::iterator iter =
        FieldFactory::data.find( solverType );

    if ( iter == FieldFactory::data.end() )
    {
        FieldFactory::data[ solverType ] =
            std::make_unique< FieldManager >();
    }
}

FieldManager * FieldFactory::GetFieldManager( int solverType )
{
    std::map< int, std::unique_ptr< FieldManager > >::iterator iter =
        FieldFactory::data.find( solverType );

    if ( iter == FieldFactory::data.end() )
    {
        return nullptr;
    }

    return iter->second.get();
}

void FieldFactory::FreeFieldManager()
{
    FieldFactory::data.clear();
}

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFaces = interFace->nIFaces;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataSend = interFace->dataSend[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
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

void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
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

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            int iBFace = grid->interFace->i2b[ iFace ];
            int tId = grid->faceTopo->rCells[ iBFace ];

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
}


void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
}

EndNameSpace
