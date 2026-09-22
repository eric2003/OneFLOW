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

class NameValuePair
{
private:
    StringField nameList;
    RealField valueList;

public:
    void Add(
        const std::string & name,
        Real value );

    void AddName(
        const std::string & name );

    int Size() const;

    const std::string & GetName( int index ) const;
    Real GetValue( int index ) const;
};

class FieldAlloc
{
public:
    static void AllocateAllFields( int solverType, const std::string & basicString );
    static void InitField( int solverType, const std::string & basicString );
    static void RegisterInterfaceVar( int solverType, const std::string & basicString );
    static void AllocateGlobalField( int solverType, const std::string & basicString );
    static void AllocateAllKindsOfInterfaceField( int solverType );
    static void AllocateInterfaceField( IFieldProperty * iFieldProperty );
    static void AllocateOversetInterfaceField( IFieldProperty * iFieldProperty );
};

class ParaNameDim
{
private:
    StringField nameList;
    IntField nEquList;

public:
    void Add(
        const std::string & name,
        int nEqu );

    int Size() const;

    const std::string & GetName( int index ) const;
    int GetNEqu( int index ) const;

    const StringField & GetNameList() const;
};

class ParaNameDimData
{
private:
    ParaNameDim comPara;
    ParaNameDim strPara;
    ParaNameDim unsPara;

public:
    ParaNameDim * GetParaNameDim( FieldCategory category );
    const ParaNameDim * GetParaNameDim( FieldCategory category ) const;
};


class FieldManager;

class ReadSuperPara
{
private:
    ParaNameDimData paraNameDimData;
    int solverType;

    FieldManager * GetFieldManager() const;

    void AddFieldProperties(
        FieldLocation location );

    void AddUnsteadyInnerFieldProperty();

public:
    explicit ReadSuperPara( int solverType );
    ~ReadSuperPara() = default;

    void Register(
        const std::string & fileName,
        FieldLocation location,
        bool isUnsteady );
};

class TextFileParser;

class BoolIO
{
private:
    StringField boolNameList;
    BoolField boolValueList;
    NameValuePair nameValuePair;

public:
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
