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

BeginNameSpace( ONEFLOW )

// ------------------------------------------------------------
// Registration categories
//
// MESG_FUNC : high-level operation orchestration / planning
// COMM_FUNC : concrete operation execution
// RECV_FUNC : communication receive handling
// TASK_FUNC : task construction / task strategy
// FILE_FUNC : file and resource preparation
// ------------------------------------------------------------

const int COMM_FUNC = 0;
const int RECV_FUNC = 1;
const int MESG_FUNC = 2;
const int TASK_FUNC = 3;
const int FILE_FUNC = 4;


class HXClone;

HXClone * GetClass(
    int operationId,
    int solverType,
    int funcType );


// ============================================================
// Operation planning
// ============================================================

void GenerateCmdList( int operationId );

// ============================================================
// Command construction
// ============================================================

void AddCmdToList( const std::string & operationName );
void AddCmdToList( int operationId, int solverType );


// ============================================================
// Task construction
// ============================================================
class Task;
Task * CreateTask( int operationId, int solverType );

// ============================================================
// Resource preparation
// ============================================================

void SetFile( Task * task, int operationId, int solverType );

// ============================================================
// Action dispatch
// ============================================================

void SetTaskAction( Task * task );

void CmdBasicAction( int funcType );
void CmdAction();
void CmdActionNext();

// ============================================================
// Operation execution entry
// ============================================================

void SingleSolverSingleGridTask( const std::string & taskName );
void MultiSolverMultiGridTask( const std::string & taskName );


EndNameSpace
