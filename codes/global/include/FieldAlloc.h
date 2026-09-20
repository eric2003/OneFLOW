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
#include <memory>
#include <string>


BeginNameSpace( ONEFLOW )

class IFieldProperty;

struct NameValuePair
{
public:
    StringField nameList;
    RealField valueList;
};

void SetFieldValues(
    int solverType,
    const NameValuePair & valuePair );

class FieldAlloc
{
public:
    static void AllocateAllFields( int solverType, const std::string & basicString );
    static void InitField( int solverType, const std::string & basicString );
    static void RegisterInterfaceVar( int solverType, const std::string & basicString );
    static void AllocateGlobalField( int solverType, const std::string & basicString );
    static void CalcInnerFieldFileName( const std::string & basicString, StringField & fileNameList );
    static void CalcInterfaceFileName( const std::string & basicString, StringField & fileNameList );
    static void CalcInterfaceFileType( IntField & fieldTypeList );
    static void AllocateAllKindsOfInterfaceField( int solverType );
    static void AllocateInterfaceField( IFieldProperty * iFieldProperty );
    static void AllocateOversetInterfaceField( IFieldProperty * iFieldProperty );
};

void AddInterfaceFieldNames(
    int solverType,
    int fieldType,
    const StringField & nameList );

struct ParaNameDim
{
public:
    StringField nameList;
    IntField dimList;
};

class ParaNameDimData
{
public:
    ParaNameDim comPara;
    ParaNameDim strPara;
    ParaNameDim unsPara;

    ParaNameDim * GetParaNameDim(
        const std::string & typeName );
};

class ReadSuperPara
{
public:
    ReadSuperPara();
    ~ReadSuperPara();
public:
    std::unique_ptr<ParaNameDimData> paraNameDimData;
    int solverType;
public:
    void Register( const std::string & fileName, int index );
    void AddFieldProperties( FieldLocation location );
    void AddUnsteadyInnerFieldProperty();
public:
    void AddBasicFieldProperty(
        ParaNameDim * paraNameDim,
        FieldLocation location,
        FieldCategory category );
};


class TextFileParser;

class BoolIO
{
public:
    BoolIO();
    ~BoolIO();
public:
    StringField boolNameList;
    BoolField boolValueList;
    NameValuePair nameValuePair;
public:
    void Add( const std::string & name, bool value );
    void ReadBool( TextFileParser & textFileParser );
    void ReadSuperBool( TextFileParser & textFileParser );
    bool CalcVarValue( const std::string & varName );
    void Read(
        TextFileParser & textFileParser,
        int valueFlag,
        ParaNameDimData * paraNameDimData = nullptr );
    void ReadFile(
        const std::string & fileName,
        int valueFlag = 0,
        ParaNameDimData * paraNameDimData = nullptr );
};

bool CalcBoolExp( bool var1, const std::string & opName, bool var2 );
bool CalcBoolExp( const std::string & varName1, const std::string & opName, const std::string & varName2 );
bool CalcVarValue( const std::string & varName, StringField & boolName, BoolField & boolVar );
int GetVarDimension( const std::string & dimName );

EndNameSpace
