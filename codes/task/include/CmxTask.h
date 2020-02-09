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


#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

const int COMM_FUNC = 0;
const int RECV_FUNC = 1;
const int MESG_FUNC = 2;
const int TASK_FUNC = 3;
const int FILE_FUNC = 4;

void CmdBasicAction( int funcType );
void CmdAction();
void CmdActionNext();

void SetTaskAction();

void GenerateCmdList( int msgId );

void AddCmdToList( const string & msgName );
void AddCmdToList( int taskCode, int solverCode );

class HXClone;
HXClone * GetClass( int msgId, int sTid, int msgType );

void CreateTask( int msgId, int sTid );
void SetFile( int msgId, int sTid );

void SsSgTask( const string & taskName );
void MsMgTask( const string & taskName );


EndNameSpace