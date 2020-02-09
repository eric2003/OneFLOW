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

#include "InterField.h"
#include "Grid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "InterFace.h"
#include "FieldRecord.h"
#include "RegisterUtil.h"
#include "DataBase.h"
#include "ActionState.h"
#include "FieldImp.h"

BeginNameSpace( ONEFLOW )

void PrepareInterfaceFieldRecord( int sTid, int iFk, int iSr, FieldRecord * fieldRecord )
{
    Grid * grid = Zone::GetGrid();
    InterFace * interFace = grid->interFace;

    InterFaceState::interFace = interFace;

    HXVector< DataStorage * > * iDataStorageList = new HXVector< DataStorage * >;

    GetInterfaceDataStorageList( iDataStorageList, iSr );

    VarNameSolver * varNameSolver = VarNameFactory::GetVarNameSolver( sTid, iFk );

    for ( int dataId = 0; dataId < iDataStorageList->size(); ++ dataId )
    {
        DataStorage * dataStorage = ( * iDataStorageList )[ dataId ];
        AddFieldRecord( fieldRecord, dataStorage, varNameSolver->data );
    }

    delete iDataStorageList;
}

void GetInterfaceDataStorageList( HXVector< DataStorage * > * iDataStorageList, int srFlag )
{
    InterFace * interFace = InterFaceState::interFace;
    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataStorage = GetInterfaceDataStorage( interFace, srFlag, ghostId );
        iDataStorageList->push_back( dataStorage );
    }
}

DataStorage * GetInterfaceDataStorage( InterFace * interFace, int srFlag, int ghostId )
{
    if ( srFlag == SEND_STORAGE )
    {
        return interFace->dataSend[ ghostId ];
    }
    else if ( srFlag == RECV_STORAGE )
    {
        return interFace->dataRecv[ ghostId ];
    }
    else
    {
        return 0;
    }
}

void AddFieldRecord( FieldRecord * fieldRecord, DataStorage * dataStorage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, filedName );
        int nEqu = GFieldProperty::GetNEqu( filedName );
        fieldRecord->AddField( field , nEqu );
    }
}

void SetInterfaceFieldData( int iSr, FieldRecord * fieldRecord )
{
    Grid * grid = Zone::GetGrid();
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int oppoSr = GetOppositeSendRecv( iSr );

    //按照设计当前zone是第zone i的第j个邻居。
    //在这里需找出zone i是当前zone的第几个邻居？这个值就是neiId。

    int neiId = interFace->z2n[ ZoneState::GetZid( oppoSr ) ];
    int nIFace = interFace->interFacePairs[ neiId ]->nIFace;
    IntField & interfaceId = interFace->GetInterfaceId( neiId, iSr );
    
    ActionState::dataBook->MoveToBegin();

    int nRecords = fieldRecord->nEquList.size();

    for ( int fieldId = 0; fieldId < nRecords; ++ fieldId )
    {
        int nEqu = fieldRecord->nEquList[ fieldId ];
        MRField * field  = fieldRecord->GetField( fieldId );
        if ( iSr == GREAT_SEND )
        {
            HXWriteSubData( ActionState::dataBook, field, interfaceId );
        }
        else
        {
            HXReadSubData( ActionState::dataBook, field, interfaceId );
        }
        
    }
}

void HXWriteSubData( DataBook * dataBook, MRField * field2D, IntField & idMap )
{
    int nEqu = field2D->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXWriteSubData( dataBook, ( * field2D )[ iEqu ], idMap );
    }
}

void HXWriteSubData( DataBook * dataBook, RealField & field, IntField & idMap )
{
    int nElem = idMap.size();
    if ( nElem <= 0 ) return;
    RealField swapField( nElem );
    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int id = idMap[ iElem ];
        swapField[ iElem ] = field[ id ];
    }
    HXWrite( dataBook, swapField );
}

void HXReadSubData( DataBook * dataBook, MRField * field2D, IntField & idMap )
{
    int nEqu = field2D->GetNEqu();
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        HXReadSubData( dataBook, ( * field2D )[ iEqu ], idMap );
    }
}

void HXReadSubData( DataBook * dataBook, RealField & field, IntField & idMap )
{
    int nElem = idMap.size();
    if ( nElem <= 0 ) return;
    RealField swapField( nElem );
    HXRead( dataBook, swapField );

    for ( int iElem = 0; iElem < nElem; ++ iElem )
    {
        int id = idMap[ iElem ];
        field[ id ] = swapField[ iElem ];
    }
}


EndNameSpace