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
#include "HXArray.h"

BeginNameSpace( ONEFLOW )

class DataBook;
class DataStorage;
class InterFace;
class FieldRecord;

DataStorage * GetInterfaceDataStorage( InterFace * interFace, int srFlag, int ghostId );
void GetInterfaceDataStorageList( HXVector< DataStorage * > * iDataStorageList, int srFlag );
void AddFieldRecord( FieldRecord * fieldRecord, DataStorage * dataStorage, StringField & fieldNameList );
void PrepareInterfaceFieldRecord( int sTid, int iFk, int iSr, FieldRecord * fieldRecord );
void SetInterfaceFieldData( int iSr, FieldRecord * fieldRecord );

void HXWriteSubData( DataBook * dataBook, MRField * field2D, IntField & idMap );
void HXWriteSubData( DataBook * dataBook, RealField & field, IntField & idMap );
void HXReadSubData( DataBook * dataBook, MRField * field2D, IntField & idMap );
void HXReadSubData( DataBook * dataBook, RealField & field, IntField & idMap );

EndNameSpace