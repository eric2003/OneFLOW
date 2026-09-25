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
#include "FieldAllocConfig.h"
#include "Fatal.h"
#include "DataBase.h"
#include "TextFileParser.h"

BeginNameSpace( ONEFLOW )

// Section A (alloc text parsing) is in FieldAllocConfig.*
namespace
{
    bool CompareValues(
        const std::string & leftToken,
        const std::string & operatorName,
        const std::string & rightToken )
    {
        int leftValue =
            ResolveIntegerValue( leftToken );

        int rightValue =
            ResolveIntegerValue( rightToken );

        if ( operatorName == ">" )
        {
            return leftValue > rightValue;
        }

        if ( operatorName == ">=" )
        {
            return leftValue >= rightValue;
        }

        if ( operatorName == "==" )
        {
            return leftValue == rightValue;
        }

        if ( operatorName == "<" )
        {
            return leftValue < rightValue;
        }

        if ( operatorName == "<=" )
        {
            return leftValue <= rightValue;
        }

        if ( operatorName == "!=" )
        {
            return leftValue != rightValue;
        }

        Fatal(
            "Unknown comparison operator: "
            + operatorName );

        return false;
    }

    bool CalcBoolLogic(
        bool leftValue,
        const std::string & operatorName,
        bool rightValue )
    {
        if ( operatorName == "&&" )
        {
            return leftValue && rightValue;
        }

        if ( operatorName == "||" )
        {
            return leftValue || rightValue;
        }

        Fatal(
            "Unknown boolean operator: "
            + operatorName );

        return false;
    }

    using BoolLineReader =
        void ( FieldConfigReader::* )( TextFileParser & );

    void ReadBoolFile(
        FieldConfigReader & configReader,
        const std::string & fileName,
        BoolLineReader trueReader )
    {
        // \t is the tab key
        std::string separator = " \r\n\t#$,;\"()";

        TextFileParser textFileParser;

        textFileParser.OpenFile(
            fileName,
            std::ios_base::in );

        textFileParser.SetDefaultSeparator(
            separator );

        while ( ! textFileParser.ReachTheEndOfFile() )
        {
            bool flag =
                textFileParser.ReadNextNonEmptyLine();

            if ( ! flag ) break;

            std::string keyWord =
                textFileParser.ReadNextWord();

            if ( keyWord == "true" )
            {
                ( configReader.*trueReader )(
                    textFileParser );
            }
            else if ( keyWord == "bool" )
            {
                configReader.ReadBool(
                    textFileParser );
            }
            else if ( keyWord == "superbool" )
            {
                configReader.ReadSuperBool(
                    textFileParser );
            }
            else
            {
                bool flag =
                    configReader.GetBoolValue( keyWord );

                if ( flag )
                {
                    ( configReader.*trueReader )(
                        textFileParser );
                }
            }
        }

        textFileParser.CloseFile();
    }
}

int ResolveIntegerValue(
    const std::string & valueToken )
{
    if ( Word::IsDigit( valueToken ) )
    {
        return StringToDigit< int >( valueToken );
    }

    return GetDataValue< int >( valueToken );
}

void FieldNameList::Add(
    const std::string & name )
{
    nameList.push_back( name );
}

int FieldNameList::Size() const
{
    return nameList.size();
}

const std::string & FieldNameList::GetName(
    int index ) const
{
    return nameList[ index ];
}

void NameValuePair::Add(
    const std::string & name,
    Real value )
{
    nameList.push_back( name );
    valueList.push_back( value );
}

int NameValuePair::Size() const
{
    return nameList.size();
}

const std::string & NameValuePair::GetName(
    int index ) const
{
    return nameList[ index ];
}

Real NameValuePair::GetValue(
    int index ) const
{
    return valueList[ index ];
}

void FieldConfigReader::Add( const std::string & name, bool value )
{
    boolNameList.push_back( name );
    boolValueList.push_back( value );
}

bool FieldConfigReader::GetBoolValue(
    const std::string & varName ) const
{
    for ( int i = 0; i < boolNameList.size(); ++ i )
    {
        if ( varName == boolNameList[ i ] )
        {
            return boolValueList[ i ];
        }
    }

    Fatal( "Unknown boolean variable: " + varName );

    return false;
}

void FieldConfigReader::ReadBool( TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    textFileParser.ReadNextWord();

    std::string var1 =
        textFileParser.ReadNextWord();

    std::string opName =
        textFileParser.ReadNextWord();

    std::string var2 =
        textFileParser.ReadNextWord();

    bool boolValue =
        CompareValues(
            var1,
            opName,
            var2 );

    this->Add( varName, boolValue );
}

void FieldConfigReader::ReadSuperBool( TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    textFileParser.ReadNextWord();

    std::string var1 =
        textFileParser.ReadNextWord();

    std::string opName =
        textFileParser.ReadNextWord();

    std::string var2 =
        textFileParser.ReadNextWord();

    bool varValue1 =
        this->GetBoolValue( var1 );

    bool varValue2 =
        this->GetBoolValue( var2 );

    bool boolValue =
        CalcBoolLogic(
            varValue1,
            opName,
            varValue2 );

    this->Add( varName, boolValue );
}

void FieldConfigReader::ReadName(
    TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    fieldNameList.Add( varName );
}

const FieldNameList & FieldConfigReader::GetFieldNameList() const
{
    return fieldNameList;
}

const NameValuePair & FieldConfigReader::GetNameValuePair() const
{
    return nameValuePair;
}

void FieldConfigReader::ReadNameValue(
    TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    Real varValue =
        textFileParser.ReadNextDigit< Real >();

    nameValuePair.Add(
        varName,
        varValue );
}

void FieldConfigReader::ReadFile(
    const std::string & fileName )
{
    ReadBoolFile(
        *this,
        fileName,
        &FieldConfigReader::ReadName );
}

void FieldConfigReader::ReadValueFile(
    const std::string & fileName )
{
    ReadBoolFile(
        *this,
        fileName,
        &FieldConfigReader::ReadNameValue );
}

EndNameSpace