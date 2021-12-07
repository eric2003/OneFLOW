/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

class IFieldProperty;

class NameValuePair
{
public:
    StringField nameList;
    RealField valueList;
};

class FieldNamePair
{
public:
    FieldNamePair();
    ~FieldNamePair();
public:
    static void SetField( int sTid, NameValuePair & valuePair );
};

class FieldAlloc
{
public:
    FieldAlloc();
    ~FieldAlloc();
public:
    static void AllocateAllFields( int sTid, const std::string & basicString );
    static void InitField( int sTid, const std::string & basicString );
    static void RegisterInterfaceVar( int sTid, const std::string & basicString );
    static void AllocateGlobalField( int sTid, const std::string & basicString );
    static void CalcInnerFieldFileName( const std::string & basicString, StringField & fileNameList );
    static void CalcInterfaceFileName( const std::string & basicString, StringField & fileNameList );
    static void CalcInterfaceFileType( IntField & fieldTypeList );
    static void AllocateAllKindsOfInterfaceField( int sTid );
    static void AllocateInterfaceField( IFieldProperty * iFieldProperty );
    static void AllocateOversetInterfaceField( IFieldProperty * iFieldProperty );
};

class ReadInterfaceVar
{
public:
    ReadInterfaceVar();
    ~ReadInterfaceVar();
public:
    static void AddFieldName( int sTid, int fieldType, StringField & nameList );
};

class ParaNameDim
{
public:
    ParaNameDim();
    ~ParaNameDim();
public:
    StringField nameList;
    IntField dimList;
};

class ParaNameDimData
{
public:
    ParaNameDimData();
    ~ParaNameDimData();
public:
    ParaNameDim * comPara;
    ParaNameDim * strPara;
    ParaNameDim * unsPara;
public:
    ParaNameDim * GetParaNameDim( const std::string & typeName );
};

class ReadSuperPara
{
public:
    ReadSuperPara();
    ~ReadSuperPara();
public:
    ParaNameDimData * paraNameDimData;
    int sTid;
public:
    void Register( const std::string & fileName, int index );
    void AddUnsteadyInnerFieldProperty();
    void AddInnerFieldProperty();
    void AddFaceFieldProperty();
    void AddBoundaryFieldProperty();
public:
    void AddBasicFieldProperty( ParaNameDim * paraNameDim, int fieldType, int type );
};


class FileIO;
class BoolIO
{
public:
    BoolIO();
    ~BoolIO();
public:
    StringField boolNameList;
    BoolField boolValueList;
    int valueFlag;
    FileIO * ioFile;
    NameValuePair nameValuePair;
    ParaNameDimData * paraNameDimData;
public:
    void Add( const std::string & name, bool value );
    void ReadBool( FileIO * ioFile );
    void ReadSuperBool( FileIO * ioFile );
    bool CalcVarValue( const std::string & varName );
    void Read();
    void ReadFile( const std::string & fileName, int valueFlag = 0 );

};

bool CalcBoolExp( bool var1, const std::string & opName, bool var2 );
bool CalcBoolExp( const std::string & varName1, const std::string & opName, const std::string & varName2 );
bool CalcVarValue( const std::string & varName, StringField & boolName, BoolField & boolVar );
int GetVarDimension( const std::string & dimName );

EndNameSpace
