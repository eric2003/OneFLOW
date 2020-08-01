/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "FileIO.h"
#include "FileUtil.h"
#include "Prj.h"

BeginNameSpace( ONEFLOW )

BcData bcdata;

BcData::BcData()
{
}

BcData::~BcData()
{
}

void BcData::Init()
{
    this->ReadList();
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
    fstream file;
    string fileName = "grid/bcRegionMap.txt";
    OpenPrjFile( file, fileName, ios_base::in );

    file >> nRegion;

    CloseFile( file );
}

void BcData::ReadList()
{
    //\tÎªtab¼ü
    string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;
    string fileName = "script/bc.txt";
    ioFile.OpenPrjFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        int regionId = ioFile.ReadNextDigit< int >();
        int num = ioFile.ReadNextDigit< int >();
        RealField f;
        this->irList.push_back( regionId );
        for ( int i = 0; i < num; ++ i )
        {
            Real value = ioFile.ReadNextDigit< Real >();
            f.push_back( value );
        }
        this->dataList.push_back( f );
    }

    ioFile.CloseFile();
}

EndNameSpace