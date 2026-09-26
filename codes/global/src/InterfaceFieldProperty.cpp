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
#include "InterfaceFieldProperty.h"
#include "FieldBase.h"
#include "Fatal.h"
#include "Grid.h"
#include "GridState.h"
#include "DataBase.h"
#include "DataStorage.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "InterFace.h"
#include "FaceTopo.h"

BeginNameSpace( ONEFLOW )

void InterfaceFieldProperty::AddField(
    const std::string & fieldName,
    int nEqu )
{
    this->fieldDefinitions.AddField(
        fieldName,
        nEqu );
}

bool InterfaceFieldProperty::HasField(
    const std::string & fieldName ) const
{
    return this->fieldDefinitions.HasField(
        fieldName );
}

int InterfaceFieldProperty::GetNEqu(
    const std::string & fieldName ) const
{
    return this->fieldDefinitions.GetNEqu(
        fieldName );
}

bool InterfaceFieldProperty::Empty() const
{
    return this->fieldDefinitions.Empty();
}

const FieldDefinitionTable::Data &
InterfaceFieldProperty::GetData() const
{
    return this->fieldDefinitions.GetData();
}

void InterfaceFieldProperty::Dump(
    std::ostream & output ) const
{
    this->fieldDefinitions.Dump( output );
}

void InterfaceFieldProperty::UploadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    const auto & data = this->GetData();
    for ( const auto & [ fieldName, nEqu ] : data )
    {
        MRField * targetField =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( targetField == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface upload: "
                + fieldName );
        }

        ONEFLOW::UploadInterfaceValue(
            grid,
            targetField,
            fieldName,
            nEqu );
    }
}

void InterfaceFieldProperty::DownloadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    const auto & data = this->GetData();
    for ( const auto & [ fieldName, nEqu ] : data )
    {
        MRField * targetField =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( targetField == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface download: "
                + fieldName );
        }

        ONEFLOW::DownloadInterfaceValue(
            grid,
            targetField,
            fieldName,
            nEqu );
    }
}

void InterfaceFieldProperty::UploadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        const FieldDefinitionTable::Data & data = this->GetData();
        for ( FieldDefinitionTable::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::UploadOversetValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void InterfaceFieldProperty::DownloadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        const FieldDefinitionTable::Data & data = this->GetData();
        for ( FieldDefinitionTable::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::DownloadOversetValue( grid, targetField, iter->first, nEqu );
        }
    }
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

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface send field is not allocated: "
                + name );
        }

        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetSId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * fieldStorage )[ iEqu ][ iFace ] =
                    ( * field2D )[ iEqu ][ iCell ];
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

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface recv field is not allocated: "
                + name );
        }

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] =
                    ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    // Reserved: overset interface transfer is not implemented yet.
    (void) grid;
    (void) field2D;
    (void) name;
    (void) nEqu;
}

void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    // Reserved: overset interface transfer is not implemented yet.
    (void) grid;
    (void) field2D;
    (void) name;
    (void) nEqu;
}

EndNameSpace