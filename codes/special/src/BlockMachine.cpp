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
    string word = ioFile->ReadNextWord();
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
    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        blkFaceSolver.DumpBlkScript();
        blkFaceSolver.SetBoundary();
        blkFaceSolver.BuildBlkFace();
        blkFaceSolver.ConstructBlockInfo();
        blkFaceSolver.DumpBcInp();
        blkFaceSolver.GenerateLineMesh();
        blkFaceSolver.GenerateFaceMesh();
        blkFaceSolver.GenerateBlkMesh();
        blkFaceSolver.DumpStandardGrid();
    }
    else
    {
        blkFaceSolver.SetBoundary();
        blkFaceSolver.BuildBlkFace();
        blkFaceSolver.ConstructBlockInfo();
        blkFaceSolver.DumpBcInp();
        blkFaceSolver.GenerateLineMesh();
        blkFaceSolver.GenerateFaceMesh();
        blkFaceSolver.GenerateBlkMesh();
        blkFaceSolver.DumpStandardGrid();

    }
}

void BlockMachine::ConstructBlockTopo()
{
}

void BlockMachine::GenerateAllBlockMesh()
{
}

void BlockMachine::DumpStandardGrid()
{
}

void BlockMachine::DumpStandardGrid( Grids & strGridList )
{
    fstream file;
    OpenPrjFile( file, "grid/strplot3d.grd", ios_base::out | ios_base::binary );

    int nZone = strGridList.size();
    HXWrite( & file, nZone );
    for ( int iBlock = 0; iBlock < nZone; ++ iBlock )
    {
        Grid * gridstr = strGridList[ iBlock ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;
        HXWrite( & file, ni );
        HXWrite( & file, nj );
        HXWrite( & file, nk );
    }

    for ( int iBlock = 0; iBlock < nZone; ++ iBlock )
    {
        Grid * gridstr = strGridList[ iBlock ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        HXWrite( & file, grid->nodeMesh->xN );
        HXWrite( & file, grid->nodeMesh->yN );
        HXWrite( & file, grid->nodeMesh->zN );
    }

    CloseFile( file );
}

void BlockMachine::DumpStandardGridBc( Grids & strGridList )
{
    fstream file;
    OpenPrjFile( file, "grid/strplot3d.inp", ios_base::out );

    int nZone = strGridList.size();
    int flowSolverIndex = 1;

    int width = 5;

    file << setw( width ) << flowSolverIndex << endl;
    file << setw( width ) << nZone << endl;

    for ( int iBlock = 0; iBlock < nZone; ++ iBlock )
    {
        Grid * gridstr = strGridList[ iBlock ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        file << setw( width ) << ni;
        file << setw( width ) << nj;
        if ( Dim::dimension == ONEFLOW::THREE_D )
        {
            file << setw( width ) << nk;
        }
        file << endl;

        string blockName = "zone";

        BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;
        
        HXVector< BcRegion * > & regions = * bcRegionGroup->regions;

        int nBcRegions = regions.size();

        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            int imin, imax, jmin, jmax, kmin, kmax;

            BcRegion * region = regions[ ir ];
            imin = region->s->start[ 0 ];
            jmin = region->s->start[ 1 ];
            kmin = region->s->start[ 2 ];

            imax = region->s->end[ 0 ];
            jmax = region->s->end[ 1 ];
            kmax = region->s->end[ 2 ];

            int bcType = region->bcType;

            file << setw( width ) << imin;
            file << setw( width ) << imax;
            file << setw( width ) << jmin;
            file << setw( width ) << jmax;
            if ( Dim::dimension == ONEFLOW::THREE_D )
            {
                file << setw( width ) << kmin;
                file << setw( width ) << kmax;
            }
            file << setw( width ) << bcType;
            file << "\n";

            if ( bcType < 0 )
            {
                imin = region->t->start[ 0 ];
                jmin = region->t->start[ 1 ];
                kmin = region->t->start[ 2 ];

                imax = region->t->end[ 0 ];
                jmax = region->t->end[ 1 ];
                kmax = region->t->end[ 2 ];

                int zid = region->t->zid;

                file << setw( width ) << imin;
                file << setw( width ) << imax;
                file << setw( width ) << jmin;
                file << setw( width ) << jmax;
                if ( Dim::dimension == ONEFLOW::THREE_D )
                {
                    file << setw( width ) << kmin;
                    file << setw( width ) << kmax;
                }
                file << setw( width ) << bcType;
                file << "\n";
            }
        }
    }

    CloseFile( file );
}

EndNameSpace