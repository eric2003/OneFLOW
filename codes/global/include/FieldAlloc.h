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

struct NameValuePair
{
public:
    StringField nameList;
    RealField valueList;
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

struct ParaNameDim
{
public:
    StringField nameList;
    IntField nEquList;
};

class ParaNameDimData
{
public:
    ParaNameDim * GetParaNameDim( FieldCategory category );
    const ParaNameDim * GetParaNameDim( FieldCategory category ) const;

public:
    ParaNameDim comPara;
    ParaNameDim strPara;
    ParaNameDim unsPara;
};

class ReadSuperPara
{
public:
    ReadSuperPara() = default;
    ~ReadSuperPara() = default;

public:
    ParaNameDimData paraNameDimData;
    int solverType;

public:
    void Register(
        const std::string & fileName,
        FieldLocation location,
        bool isUnsteady );

    void AddFieldProperties( FieldLocation location );
    void AddUnsteadyInnerFieldProperty();
};

class TextFileParser;

class BoolIO
{
private:
    StringField boolNameList;
    BoolField boolValueList;

public:
    NameValuePair nameValuePair;

public:
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
