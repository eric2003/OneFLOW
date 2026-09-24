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
#include "HXDefine.h"
#include "FieldCategory.h"
#include <string>


BeginNameSpace( ONEFLOW )

class IFieldProperty;
struct UsdFieldNames;

class FieldNameList
{
private:
    StringField nameList;

public:
    void Add(
        const std::string & name );

    int Size() const;

    const std::string & GetName(
        int index ) const;
};

class NameValuePair
{
private:
    StringField nameList;
    RealField valueList;

public:
    void Add(
        const std::string & name,
        Real value );

    int Size() const;

    const std::string & GetName(
        int index ) const;

    Real GetValue(
        int index ) const;
};

class FieldManager;
class UnsGrid;
class FieldPropertyData;

class FieldAlloc
{
public:
    static void AllocateAllFields(
        int solverType,
        const std::string & basicString );
};

class TextFileParser;

class BoolIO
{
private:
    StringField boolNameList;
    BoolField boolValueList;

    FieldNameList fieldNameList;
    NameValuePair nameValuePair;

public:
    const FieldNameList & GetFieldNameList() const;

    const NameValuePair & GetNameValuePair() const;

    bool GetBoolValue(
        const std::string & varName ) const;

    void Add(
        const std::string & name,
        bool value );

    void ReadBool(
        TextFileParser & textFileParser );

    void ReadSuperBool(
        TextFileParser & textFileParser );

    void ReadName(
        TextFileParser & textFileParser );

    void ReadNameValue(
        TextFileParser & textFileParser );

    void ReadFile(
        const std::string & fileName );

    void ReadValueFile(
        const std::string & fileName );
};

EndNameSpace
