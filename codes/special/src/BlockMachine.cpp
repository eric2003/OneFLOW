/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "BlockMachine.h"
#include "LineMachine.h"
#include "CurveInfo.h"
#include "FileIO.h"
#include "FileUtil.h"
#include "BgGrid.h"
#include "StrGrid.h"
#include "GridState.h"
#include "NodeMesh.h"
#include "HXPointer.h"
#include "DataBaseIO.h"
#include "Prj.h"
#include "BcRecord.h"
#include "Dimension.h"
#include "BlockElem.h"
#include "BlockFaceSolver.h"
#include "GridPara.h"
#include "HXCgns.h"
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

BlockMachine block_Machine;

BlockMachine::BlockMachine()
{
}

BlockMachine::~BlockMachine()
{
}

void BlockMachine::AddFaceToBlock( FileIO * ioFile )
{
    std::string word = ioFile->ReadNextWord();
    if ( word == "L2F" )
    {
        int faceid = ioFile->ReadNextDigit< int >();
        int pos = ioFile->ReadNextDigit< int >();
        int lineid = ioFile->ReadNextDigit< int >();
        
        if ( Dim::dimension == ONEFLOW::THREE_D )
        {
            blkFaceSolver.AddLineToFace( faceid, pos, lineid );
        }
        else
        {
            blkFaceSolver.AddLineToFace( faceid, pos, lineid );
        }
    }
    else if ( word == "F2B" )
    {
        int blockid = ioFile->ReadNextDigit< int >();
        int pos = ioFile->ReadNextDigit< int >();
        int faceid = ioFile->ReadNextDigit< int >();

        if ( Dim::dimension == ONEFLOW::THREE_D )
        {
            blkFaceSolver.AddFace2Block( blockid, pos, faceid );
        }
        else
        {
            blkFaceSolver.AddFace2Block( blockid, pos, faceid );
        }
    }
}

void BlockMachine::GenerateFaceBlockLink()
{
    blkFaceSolver.GenerateFaceBlockLink();
}

EndNameSpace
