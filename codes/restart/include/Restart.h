/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class Restart
{
public:
    Restart();
    virtual ~Restart();
public:
    int sTid;
public:
    void ReadUnsteady( int sTid );
    void DumpUnsteady( int sTid );
    void InitUnsteady( int sTid );
    void Read( int sTid );
	//void Readins(int sTid);
    void Dump( int sTid );
public:
    virtual void InitRestart( int sTid );
	virtual void InitinsRestart( int sTid );
};

Restart * CreateRestart( int sTid );

class DataStorage;

void ReadRestartHeader();
void ReadinsRestartHeader();
void DumpRestartHeader();
void RwInterface( int sTid, int readOrWrite );
void RwInterfaceRecord( DataStorage * storage, StringField & fieldNameList, int readOrWrite );
void ReadFieldRecord( DataStorage * storage, StringField & fieldNameList );
void WriteFieldRecord( DataStorage * storage, StringField & fieldNameList );
void ZeroFieldRecord( DataStorage * storage, StringField & fieldNameList );

EndNameSpace