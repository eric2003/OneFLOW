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
#include "GridDef.h"

BeginNameSpace( ONEFLOW )

class BlockInfo;
class BlockMesh;
class Grid;
class FileIO;

class BlockMachine
{
public:
    BlockMachine();
    ~BlockMachine();
public:
    HXVector< BlockInfo * > blockInfoList;
    HXVector< BlockMesh * > blockMeshList;
public:
    void AddFaceToBlock( FileIO * ioFile );
    void ConstructBlockTopo();
    void GenerateAllBlockMesh();
    void DumpStandardGrid();
    void DumpStandardGrid( Grids & strGridList );
    void DumpStandardGridBc( Grids & strGridList );
    void GenerateFaceBlockLink();
};

extern BlockMachine block_Machine;


EndNameSpace