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

#include "BcData.h"
#include "TextFileParser.h"

#include "Prj.h"

BeginNameSpace( ONEFLOW )

BcData::BcData()
{
}

BcData::~BcData()
{
}

void BcData::Init( const std::string & fileName )
{
    this->ReadList( fileName );
    this->ReadRegion();
    this->r2d.resize( nRegion, -1 );
    for ( int i = 0; i < irList.size(); ++ i )
    {
        int ir = irList[ i ];
        this->r2d[ ir ] = i;
    }
}

void BcData::ReadRegion()
{
    std::fstream file;
    std::string fileName = "grid/bcRegionMap.txt";
    Prj::OpenPrjFile( file, fileName, std::ios_base::in );

    file >> nRegion;

    Prj::CloseFile( file );
}

void BcData::ReadList( const std::string & fileName )
{
    //\t is the tab key
    std::string separator = " =\r\n\t#$,;\"";

    TextFileParser textFileParser;
    textFileParser.OpenPrjFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        int regionId = textFileParser.ReadNextDigit< int >();
        int num = textFileParser.ReadNextDigit< int >();
        RealField f;
        this->irList.push_back( regionId );
        for ( int i = 0; i < num; ++ i )
        {
            Real value = textFileParser.ReadNextDigit< Real >();
            f.push_back( value );
        }
        this->dataList.push_back( f );
    }

    textFileParser.CloseFile();
}

EndNameSpace
