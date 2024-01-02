/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "Partition.h"
#include "metis.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "Stop.h"
#include "BcRecord.h"
#include "FaceTopo.h"
#include "CellTopo.h"
#include "CellMesh.h"
#include "BgGrid.h"
#include "GridState.h"
#include "NodeMesh.h"
#include "InterFace.h"
#include "CalcGrid.h"
#include "DataBase.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

L2GMapping::L2GMapping()
{
    ;
}

L2GMapping::~L2GMapping()
{
    ;
}

void L2GMapping::CalcL2G( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    this->Alloc( grid );
    this->CalcL2GNode( ggrid, zid, grid );
    this->CalcL2GFace( ggrid, zid, grid );
    this->CalcL2GCell( ggrid, zid, grid );
}

void L2GMapping::Alloc( UnsGrid * grid )
{
    int nNodes = grid->nNodes;
    int nFaces = grid->nFaces;
    int nCells = grid->nCells;

    this->l2g_node.resize( nNodes );
    this->l2g_face.resize( nFaces );
    this->l2g_cell.resize( nCells );
}

void L2GMapping::CalcL2GNode( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nNodes = ggrid->nNodes;

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        if ( g2l->g2l_node[ iNode ] > - 1 )
        {
            this->l2g_node[ g2l->g2l_node[ iNode ] ] = iNode;
        }
    }
}

void L2GMapping::CalcL2GFace( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nFaces = ggrid->nFaces;

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int fid = g2l->g2l_face[ iFace ];
        if ( fid >= 0 )
        {
            this->l2g_face[ fid ] = iFace;
        }
    }
}

void L2GMapping::CalcL2GCell( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nCells = ggrid->nCells;
    int cid = 0;
    for ( int gcid = 0; gcid < nCells; ++ gcid )
    {
        if ( g2l->gc2lzone[ gcid ] == zid )
        {
            this->l2g_cell[ cid ] = gcid;
            ++ cid;
        }
    }
}

G2LMapping::G2LMapping( UnsGrid * ggrid )
{
    this->ggrid = ggrid;

    this->g2l_cell.resize( ggrid->nCells );
    this->g2l_face.resize( ggrid->nFaces );
    this->g2l_node.resize( ggrid->nNodes );
    this->gc2lzone.resize( ggrid->nCells );

    this->npartproc = GetDataValue< int >( "npartproc" );
}

G2LMapping::~G2LMapping()
{
}

void G2LMapping::GenerateGC2Z()
{
    if ( npartproc < 2 )
    {
        Stop( "The number of partitions should be greater than 1!\n" );
    }

    int nCells  = ggrid->nCells;
    int nFaces  = ggrid->nFaces;
    int nBFaces = ggrid->nBFaces;

    std::vector<idx_t> xadj  ( ggrid->nCells + 1 );
    std::vector<idx_t> adjncy( 2 * ( nFaces - nBFaces ) );

    this->GetXadjAdjncy( ggrid, xadj, adjncy );
    this->PartByMetis( nCells, xadj, adjncy );
    //this->DumpGC2Z( gridForPartition );
    //this->ReadGC2Z( gridForPartition );
}
#ifdef ENABLE_METIS
void G2LMapping::GetXadjAdjncy( UnsGrid * ggrid, std::vector<idx_t> & xadj, std::vector<idx_t>& adjncy )
{   
    int  nCells = ggrid->nCells;
    CalcC2C( ggrid );
    LinkField & c2c = ggrid->cellMesh->cellTopo->c2c;
    xadj[ 0 ]  = 0;
    int iCount = 0;
    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        xadj[ iCell + 1 ] = xadj[ iCell ] + c2c[ iCell ].size();
        for ( HXSize_t j = 0; j < c2c[ iCell ].size(); ++ j )
        {
            adjncy[ iCount ++ ] = c2c[ iCell ][ j ];
        }
    }
}

void G2LMapping::PartByMetis( idx_t nCells, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy )
{
    idx_t   ncon     = 1;
    idx_t   * vwgt   = 0;
    idx_t   * vsize  = 0;
    idx_t   * adjwgt = 0;
    float * tpwgts = 0;
    float * ubvec  = 0;
    idx_t options[ METIS_NOPTIONS ];
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t objval;
    idx_t nZone = npartproc;

    METIS_SetDefaultOptions( options );
    std::cout << "Now begining partition graph!\n";
    if ( nZone > 8 )
    {
        std::cout << "Using K-way Partitioning!\n";
        METIS_PartGraphKway( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
                             & nZone, tpwgts, ubvec, options, & objval, & gc2lzone[ 0 ] );
    }
    else
    {
        std::cout << "Using Recursive Partitioning!\n";
        METIS_PartGraphRecursive( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
                                  & nZone, tpwgts, ubvec, options, & objval, & gc2lzone[ 0 ] );
    }
    std::cout << "The interface number: " << objval << std::endl; 
    std::cout << "Partition is finished!\n";
}
#endif

void G2LMapping::DumpXadjAdjncy( UnsGrid * grid, IntField & xadj, IntField & adjncy )
{
    ;
}

void G2LMapping::DumpGC2Z( UnsGrid * grid )
{
    ;
}

void G2LMapping::ReadGC2Z( UnsGrid * grid )
{
    ;
}


Partition::Partition()
{
    g2l = 0;
    this->partition_type = GetDataValue< int >( "partition_type" );
    this->npartproc = GetDataValue< int >( "npartproc" );
    this->partition_c2n = GetDataValue< int >( "partition_c2n" );
}

Partition::~Partition()
{
    delete g2l;
}

void Partition::Run()
{
    this->ReadGrid();

    this->GenerateMultiZoneGrid();

    ONEFLOW::GenerateMultiZoneCalcGrids( grids );
}

void Partition::ReadGrid()
{
    StringField gridFileList;
    std::string ori_uns_file = GetDataValue< std::string >( "ori_uns_file" );
    gridFileList.push_back( ori_uns_file );

    Zone::ReadGrid( gridFileList );

    int nZones = ZoneState::nZones;

    if ( nZones > 1 )
    {
        std::cout << " At present, there is no support for multiple blocks such as nZones > 1 !\n";
    }
    else
    {
        Grid * grid = Zone::GetGrid();
        uns_grid = UnsGridCast( grid );
    }
}

void Partition::GenerateMultiZoneGrid()
{
    this->CreatePart();

    this->AllocPart();

    this->BuildCalculationalGrid();
}

void Partition::CreatePart()
{
    g2l = new G2LMapping( uns_grid );
    g2l->npartproc = this->npartproc;
    g2l->GenerateGC2Z();
    this->CalcG2lCell();

    if ( this->partition_c2n )
    {
        this->CalcGC2N();
    }
}

void Partition::AllocPart()
{
    grids.resize( npartproc );
    for ( int pid = 0; pid < npartproc; ++ pid )
    {
        int gridType = ONEFLOW::UMESH;
        Grid * grid = ONEFLOW::CreateGrid( gridType );
        grid->level = 0;
        grid->id = pid;
        grid->localId = pid;
        grid->type = gridType;
        grids[ pid ] = grid;
    }
}

void Partition::BuildCalculationalGrid()
{
    this->PreProcess();
    for ( int pid = 0; pid < npartproc; ++ pid )
    {
        std::cout << "BuildCalculationalGrid pid = " << pid << " npartproc = " << npartproc << "\n";
        //for unstructured grid, each processor only contains one zone, so pid equal to zid
        this->BuildCalculationalGrid( pid );
    }
    this->PostProcess();
}

void Partition::PreProcess()
{
}

void Partition::PostProcess()
{
}

void Partition::CalcGC2N()
{
}

void Partition::CalcG2lCell()
{
    IntField zCount( npartproc, 0 );

    int nCells = uns_grid->nCells;
    for ( int cid = 0; cid < nCells; ++ cid )
    {
        int zid = g2l->gc2lzone[ cid ];
        g2l->g2l_cell[ cid ] = zCount[ zid ] ++;
    }
}

void Partition::BuildCalculationalGrid( int zid )
{
    UnsGrid * grid = UnsGridCast( grids[ zid ] );

    grid->nCells = this->GetNCell( uns_grid, zid );

    this->CalcG2lFace( uns_grid, zid, grid );
    this->CalcG2lNode( uns_grid, zid, grid );

    if ( partition_c2n )
    {
    //    this->CalcGlobalToLocalCellToNodeMapping( uns_grid, zid, grid );
    //    this->WriteCellToNode( grid );
    }

    this->l2g = new L2GMapping();
    this->CreateL2g( uns_grid, zid, grid );
    this->SetCoor  ( uns_grid, zid, grid );
    this->SetGeometricRelationship( uns_grid, zid, grid );
    delete this->l2g;
}

void Partition::CalcG2lFace( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nCells  = ggrid->nCells;
    int nFaces  = ggrid->nFaces;
    int nBFaces = ggrid->nBFaces;

    IntField & glCell = ggrid->faceTopo->lCells;
    IntField & grCell = ggrid->faceTopo->rCells;

    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        g2l->g2l_face[ fid ] = - 2;
    }

    //set all face in iZone to -1
    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        int glc = glCell[ fid ];
        int grc = grCell[ fid ];
        if ( g2l->gc2lzone[ glc ] == zid )
        {
            g2l->g2l_face[ fid ] = - 1;
        }
        else if ( grc < nCells && g2l->gc2lzone[ grc ] == zid )
        {
            g2l->g2l_face[ fid ] = - 1;
        }
    }

    int nFaceNow = 0;

    //physical boundary
    for ( int fid = 0; fid < nBFaces; ++ fid )
    {
        if ( g2l->g2l_face[ fid ] == - 1 )
        {
            g2l->g2l_face[ fid ] = nFaceNow ++;
        }
    }

    int nIFaceNow = 0;
    for ( int fid = nBFaces; fid < nFaces; ++ fid )
    {
        if ( g2l->g2l_face[ fid ] == - 1 )
        {
            int glc = glCell[ fid ];
            int grc = grCell[ fid ];
            if ( g2l->gc2lzone[ glc ] != g2l->gc2lzone[ grc ] )
            {
                g2l->g2l_face[ fid ] = nFaceNow ++;
                nIFaceNow ++;
            }
        }
    }

    //inner boundary
    int nBFaceNow = nFaceNow;
    for ( int fid = nBFaces; fid < nFaces; ++ fid )
    {
        if ( g2l->g2l_face[ fid ] == - 1 )
        {
            int glc = glCell[ fid ];
            int grc = grCell[ fid ];

            if ( g2l->gc2lzone[ glc ] == g2l->gc2lzone[ grc ] )
            {
                g2l->g2l_face[ fid ] = nFaceNow ++;
            }
        }
    }

    grid->nFaces  = nFaceNow;
    grid->nBFaces = nBFaceNow;

    InterFace * interFace = grid->interFace;
    interFace->Set( nIFaceNow );
    grid->nIFaces = nIFaceNow;
}

void Partition::CalcG2lNode( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nFaces = ggrid->nFaces;
    int nNodes = ggrid->nNodes;

    LinkField & f2n = ggrid->faceTopo->faces;

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        g2l->g2l_node[ iNode ] = - 2;
    }

    //set iZone g2l->g2l_node to -1
    int iCount = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( g2l->g2l_face[ iFace ] > - 1 )
        {
            int nFNode = f2n[ iFace ].size();
            for ( int iNode = 0; iNode < nFNode; ++ iNode )
            {
                g2l->g2l_node[ f2n[ iFace ][ iNode ] ] = - 1;
            }
        }
    }

    int nLNode = 0;
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        if ( g2l->g2l_node[ iNode ] == - 1 )
        {
            g2l->g2l_node[ iNode ] = nLNode ++;
        }
    }

    grid->nNodes = nLNode;
}

int Partition::GetNCell( UnsGrid * ggrid, int zid )
{
    int nCells = ggrid->nCells;
    int iCount = 0;
    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        if ( g2l->gc2lzone[ iCell ] == zid )
        {
            iCount ++;
        }
    }
    return iCount;
}

void Partition::CreateL2g( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    l2g->g2l = this->g2l;
    l2g->CalcL2G( ggrid, zid, grid );
}

void Partition::SetCoor( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nNodes = grid->nNodes;
    grid->nodeMesh->CreateNodes( nNodes );

    int iCount = 0;
    for ( int iNode = 0; iNode < ggrid->nNodes; ++ iNode )
    {
        if ( g2l->g2l_node[ iNode ] > - 1 )
        {
            grid->nodeMesh->xN[ iCount ] = ggrid->nodeMesh->xN[ iNode ];
            grid->nodeMesh->yN[ iCount ] = ggrid->nodeMesh->yN[ iNode ];
            grid->nodeMesh->zN[ iCount ] = ggrid->nodeMesh->zN[ iNode ];
            ++ iCount;
        }
    }

    if ( iCount != nNodes )
    {
        std::cout << "error in Partition::SetCoor\n";
    }
}

void Partition::SetGeometricRelationship( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    this->CalcF2N( ggrid, zid, grid );
    this->SetF2CAndBC( ggrid, zid, grid );
    this->SetInterface( ggrid, zid, grid );
}

void Partition::CalcF2N( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    LinkField & f2n = grid->faceTopo->faces;
    LinkField & gf2n = ggrid->faceTopo->faces;

    int nFaces = grid->nFaces;
    f2n.resize( nFaces );

    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        int gfid = l2g->l2g_face[ fid ];

        int nFNode = gf2n[ gfid ].size();

        for ( int iNode = 0; iNode < nFNode; ++ iNode )
        {
            int gnid = gf2n[ gfid ][ iNode ];
            int nid  = g2l->g2l_node[ gnid ];
            f2n[ fid ].push_back( nid );
        }
    }
}

void Partition::SetF2CAndBC( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    int nGBFace = ggrid->nBFaces;

    IntField & glCell = ggrid->faceTopo->lCells;
    IntField & grCell = ggrid->faceTopo->rCells;

    IntField & gbcType = ggrid->faceTopo->bcManager->bcRecord->bcType;

    int nFaces  = grid->nFaces;
    int nBFaces = grid->nBFaces;

    IntField & lCell = grid->faceTopo->lCells;
    IntField & rCell = grid->faceTopo->rCells;
    lCell.resize( nFaces );
    rCell.resize( nFaces );

    grid->faceTopo->SetNBFaces( nBFaces );

    IntField & local_bcType = grid->faceTopo->bcManager->bcRecord->bcType;

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int gfid = this->l2g->l2g_face[ iFace ];

        int glc = glCell[ gfid ];
        int grc = grCell[ gfid ];

        int lc, rc, bctype;

        if ( gfid < nGBFace )
        {
            rc = - 1;
            lc = this->g2l->g2l_cell[ glc ];
            bctype = gbcType[ gfid ];
         }
        else
        {
            bctype = -1;
            // int face
            if ( this->g2l->gc2lzone[ glc ] == zid )
            {
                lc = this->g2l->g2l_cell[ glc ];
                rc = - 1;
            }
            else if ( this->g2l->gc2lzone[ grc ] == zid )
            {
                rc = this->g2l->g2l_cell[ grc ];
                lc = - 1;
            }
            else
            {
                std::cout << "error in SetF2CAndBC\n";
            }
        }

        local_bcType[ iFace ] = bctype;

        lCell[ iFace ] = lc;
        rCell[ iFace ] = rc;
    }

    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        int gfid = this->l2g->l2g_face[ iFace ];
        int glc = glCell[ gfid ];
        int grc = grCell[ gfid ];

        int lc = this->g2l->g2l_cell[ glc ];
        int rc = this->g2l->g2l_cell[ grc ];

        lCell[ iFace ] = lc;
        rCell[ iFace ] = rc;
    }
}

void Partition::SetInterface( UnsGrid * ggrid, int zid, UnsGrid * grid )
{
    if ( this->partition_type != 1 ) return;

    InterFace * interFace = grid->interFace;
    int nIFaces = interFace->nIFaces;
    int nBFaces = grid->nBFaces;

    int nGFace = ggrid->nFaces;
    int nGBFace = ggrid->nBFaces;

    IntField & glCell = ggrid->faceTopo->lCells;
    IntField & grCell = ggrid->faceTopo->rCells;

    //number of physical boundary face
    int nPBFace = nBFaces - nIFaces;

    for ( int gfid = nGBFace; gfid < nGFace; ++ gfid )
    {
        int fid = this->g2l->g2l_face[ gfid ];
        if ( fid < nBFaces && fid > - 1 )
        {
            //local interface id
            int ifid = fid - nPBFace;

            int glc = glCell[ gfid ];
            int grc = grCell[ gfid ];

            int leftZone  = this->g2l->gc2lzone[ glc ];
            int rightZone = this->g2l->gc2lzone[ grc ];

            int gcid = -1;

            if ( leftZone == zid )
            {
                interFace->idir[ ifid ] = 1;
                gcid = grc;
            }
            else if ( rightZone == zid )
            {
                interFace->idir[ ifid ] = - 1;
                gcid = glc;
            }
            else
            {
                std::cout << "error\n";
            }
            //interface
            //   |zone
            //   |face
            //   |cell
            //
            interFace->zoneId          [ ifid ] = this->g2l->gc2lzone[ gcid ];
            interFace->localInterfaceId[ ifid ] = this->g2l->g2l_cell[ gcid ];
            interFace->localCellId     [ ifid ] = this->g2l->g2l_cell[ gcid ];
            interFace->i2b             [ ifid ] = this->g2l->g2l_face[ gfid ];
        }
    }
}

bool FindMatch( UnsGrid * grid, FacePair * facePair )
{
    bool found = false;

    InterFace * interFace = grid->interFace;

    if ( ! ONEFLOW::IsValid( interFace ) ) return found;

    int nBFaces = grid->nBFaces;
    int nIFaces = interFace->nIFaces;
    int nPBFace = nBFaces - nIFaces;

    IntField & lCell = grid->faceTopo->lCells;
    IntField & rCell = grid->faceTopo->rCells;

    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int lc = lCell[ iFace + nPBFace ];
        int rc = rCell[ iFace + nPBFace ];
        int cell_id = MAX( lc, rc );

        if ( ( interFace->zoneId[ iFace ]      == facePair->lf.zone_id ) &&
             ( interFace->localCellId[ iFace ] == facePair->lf.cell_id ) && 
             ( cell_id                         == facePair->rf.cell_id ) )
        {
            interFace->localInterfaceId[ iFace ] = facePair->lf.face_id;
            facePair->rf.face_id = iFace;
            found = true;
            break;
        } 
    }

    return found;
}

EndNameSpace
