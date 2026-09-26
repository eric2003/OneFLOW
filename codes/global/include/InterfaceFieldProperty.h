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
#pragma once
#include "NamespaceMacros.h"
#include "FieldDefinitionTable.h"
#include "HXArray.h"
#include <string>

BeginNameSpace( ONEFLOW )

class UnsGrid;

// InterfaceFieldProperty is the Interface Storage field registry: which
// fields must exist on the interface send/recv DataStorage, and how to
// move their values to/from the grid's own fields. Unlike
// FieldDefinitionTable/FieldDefinitionSet, this class DOES touch live
// Grid / DataStorage objects (see AllocateInterfaceField and the
// Upload/Download methods below).
class InterfaceFieldProperty
{
public:
    void AddField(
        const std::string & fieldName,
        int nEqu );

    bool HasField(
        const std::string & fieldName ) const;

    int GetNEqu(
        const std::string & fieldName ) const;

    bool Empty() const;

    const FieldDefinitionTable::Data & GetData() const;

    void Dump(
        std::ostream & output ) const;
private:
    FieldDefinitionTable fieldDefinitions;
};

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );

EndNameSpace