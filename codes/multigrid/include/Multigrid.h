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

class TimeSpan;

class MG
{
public:
    typedef void ( MG::* FunctionPointer )( int );
public:
    MG();
    ~MG();
public:
    static int nPre;
    static int nPost;
    static int nMulti;
    static int iterMode;
    static int cycleType;
    static int currentGridLevel;
public:
    static void Init();
public:
    void MultigridSolve();
    void Run();
    void Allocate();
    void Deallocate();
public:
    void SolveInnerIter();
    void StrongIter();
    void WeakIter();
public:
    void MWrap( FunctionPointer multigridPointer, int gridLevel );
protected:
    void ZeroResidualsForAllSolvers();
protected:
    void PreprocessMultigridFlowField( int gl );
    void PostRelaxationCycle( int gl );
    void PreRelaxationCycle( int gl );
    void InitializeCoarseGridFlowFieldByRestrictFineGridFlowField( int fgl );
    void StoreCoarseGridFlowFieldToTemporaryStorage( int fgl );
    void PrepareFineGridResiduals( int fgl );
    void PrepareCoarseGridResiduals( int fgl );
    void CorrectFineGridFlowFieldByInterplateCoarseGridFlowField( int fgl );
    void PostprocessMultigridFlowField( int gl );
    void SolveMultigridFlowField( int gl );
    void FastSolveFlowFieldByMultigridMethod( int gl );
    void SolveCoarseGridFlowField( int fgl );
protected:
    void InnerProcess();
    void OuterProcess( TimeSpan * timeSpan );
};

bool DoNotNeedMultigridMethod( int gl );

void MultigridSolve();

EndNameSpace