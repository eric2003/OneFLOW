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
#include "BoundaryExchangeContract.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

namespace
{

// Reused endpoint scratch for the interface pack/unpack helpers.  The host
// field path is single threaded, and the ids are rebuilt per ghost level, so a
// thread-local buffer removes the per-call allocation without caching any
// topology that could go stale.
std::vector< int > & LocalInterfaceCellIds()
{
    static thread_local std::vector< int > cellIds;
    return cellIds;
}

}

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

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace.get();
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFaces = interFace->nIFaces;

    if ( field2D == nullptr ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataSend = &interFace->GetSendStorage( ghostId );

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface send field is not allocated: "
                + name );
        }

        // The send storage is already equation-major with one entry per
        // interface face, i.e. exactly the packed layout.  Copy straight into
        // it instead of packing into a temporary vector and copying again;
        // the endpoint range check is fused into the single id pass.
        const int cellCount =
            static_cast< int >( ( * field2D )[ 0 ].size() );
        std::vector< int > & cellIds = LocalInterfaceCellIds();
        cellIds.resize( static_cast< std::size_t >( nIFaces ) );
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell = -1;
            grid->GetFaceTopo().GetSId(iFace, ghostId + 1, iCell);
            if ( iCell < 0 || iCell >= cellCount )
            {
                throw std::out_of_range(
                    "boundary endpoint cell is out of range" );
            }
            cellIds[ iFace ] = iCell;
        }
        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            const Real * source = ( * field2D )[ iEqu ].data();
            Real * target = ( * fieldStorage )[ iEqu ].data();
            for ( int iFace = 0; iFace < nIFaces; ++ iFace )
            {
                target[ iFace ] = source[ cellIds[ iFace ] ];
            }
        }
    }
}

void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace.get();
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == nullptr ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = &interFace->GetRecvStorage( ghostId );

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface recv field is not allocated: "
                + name );
        }

        // Mirror of UploadInterfaceValue: the receive storage already has the
        // packed per-face layout, so scatter it straight into the ghost cells
        // of the target field without an intermediate packed copy.
        const int nIFaces = interFace->nIFaces;
        const int cellCount =
            static_cast< int >( ( * field2D )[ 0 ].size() );
        std::vector< int > & cellIds = LocalInterfaceCellIds();
        cellIds.resize( static_cast< std::size_t >( nIFaces ) );
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell = -1;
            grid->GetFaceTopo().GetTId(iFace, ghostId + 1, iCell);
            if ( iCell < 0 || iCell >= cellCount )
            {
                throw std::out_of_range(
                    "boundary endpoint cell is out of range" );
            }
            cellIds[ iFace ] = iCell;
        }
        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            const Real * source = ( * fieldStorage )[ iEqu ].data();
            Real * target = ( * field2D )[ iEqu ].data();
            for ( int iFace = 0; iFace < nIFaces; ++ iFace )
            {
                target[ cellIds[ iFace ] ] = source[ iFace ];
            }
        }
    }
}

void UploadInterfaceValuePacked(
    UnsGrid * grid, MRField * const * fields, const std::string * names,
    int fieldCount, int nEqu, int ghostId, const Real * packed )
{
    InterFace * interFace = grid->interFace.get();
    if ( ! ONEFLOW::IsValid( interFace ) ) return;
    if ( fieldCount <= 0 || fields == nullptr || names == nullptr
         || packed == nullptr )
    {
        return;
    }
    const int nIFaces = interFace->nIFaces;
    DataStorage * dataSend = &interFace->GetSendStorage( ghostId );
    for ( int fieldId = 0; fieldId < fieldCount; ++ fieldId )
    {
        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataSend, names[ fieldId ] );
        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface send field is not allocated: "
                + names[ fieldId ] );
        }
        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            const Real * source =
                packed + ( static_cast< std::size_t >( fieldId ) * nEqu + iEqu )
                    * nIFaces;
            Real * target = ( * fieldStorage )[ iEqu ].data();
            std::copy( source, source + nIFaces, target );
        }
    }
}

void DownloadInterfaceValuePacked(
    UnsGrid * grid, MRField * const * fields, const std::string * names,
    int fieldCount, int nEqu, int ghostId, Real * packed )
{
    InterFace * interFace = grid->interFace.get();
    if ( ! ONEFLOW::IsValid( interFace ) ) return;
    if ( fieldCount <= 0 || fields == nullptr || names == nullptr
         || packed == nullptr )
    {
        return;
    }
    const int nIFaces = interFace->nIFaces;
    DataStorage * dataRecv = &interFace->GetRecvStorage( ghostId );
    for ( int fieldId = 0; fieldId < fieldCount; ++ fieldId )
    {
        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataRecv, names[ fieldId ] );
        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface recv field is not allocated: "
                + names[ fieldId ] );
        }
        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            Real * target =
                packed + ( static_cast< std::size_t >( fieldId ) * nEqu + iEqu )
                    * nIFaces;
            const Real * source = ( * fieldStorage )[ iEqu ].data();
            std::copy( source, source + nIFaces, target );
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
