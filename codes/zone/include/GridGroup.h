/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "GridDef.h"
#include <fstream>
#include <string>


BeginNameSpace( ONEFLOW )

class Grid;
class UnsGrid;

class GridGroup
{
public:
    GridGroup( int zoneStart = 0 );
    ~GridGroup();
public:
    int   nZones;
    IntField pid;
    IntField zoneType;
    int zoneStart;
protected:
    void ReadGrid( std::fstream & file, int iZone );
    void CreateGrid( int zoneId );
    void CreateGridImp( int zoneId );
    void CreateGridTest( int zoneId );
public:
    void ReadGrid( const std::string & fileName );
    void InitZoneLayout( const std::string & fileName );
protected:
    void InitZoneLayout( std::fstream & file );
    void SetMultiZoneLayout();
};

class DataBook;
void ReadAbstractData( std::fstream & file, DataBook * dataBook, int sendpid, int recvpid, int tag = 0 );
void DataToGrid( DataBook * dataBook, int zid );
void DataToGridImp( DataBook * dataBook, int zid );
void DataToGridTest( DataBook * dataBook, int zid );


EndNameSpace
