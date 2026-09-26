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
#include "FieldCategory.h"
#include <map>
#include <ostream>
#include <string>

BeginNameSpace( ONEFLOW )

// FieldSpec describes the allocation definition of a field, as parsed from
// one line of a system/<solver>/alloc/*.txt file.
struct FieldSpec
{
    std::string name;
    int nEqu;
    FieldApplicability applicability;
};

// FieldDefinitionTable is a pure name -> nEqu registry. It has no dependency
// on Grid, Zone, or DataStorage: it only records what has been *defined*,
// not where it is *allocated*. Runtime allocation lives in FieldAllocator.cpp
// (Section C) and in InterfaceFieldProperty (for Interface Storage).
class FieldDefinitionTable
{
public:
    using Data = std::map< std::string, int >;

public:
    void AddField(
        const std::string & fieldName,
        int nEqu );

    bool HasField(
        const std::string & fieldName ) const;

    int GetNEqu(
        const std::string & fieldName ) const;

    bool Empty() const;

    const Data & GetData() const;

    void Dump(
        std::ostream & output ) const;

private:
    Data data;
};

// FieldDefinitionSet groups the three FieldDefinitionTable slots (Inner /
// Face / Boundary) that FieldManager keeps per FieldApplicability category.
class FieldDefinitionSet
{
public:
    FieldDefinitionTable & GetFieldDefinition(
        FieldLocation location );

    const FieldDefinitionTable & GetFieldDefinition(
        FieldLocation location ) const;

private:
    FieldDefinitionTable bcField;
    FieldDefinitionTable faceField;
    FieldDefinitionTable innerField;
};

EndNameSpace