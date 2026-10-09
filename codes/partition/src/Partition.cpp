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

#include "Partition.h"
#include <memory>
#include "metis.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "Fatal.h"
#include "BcRecord.h"
#include "FaceTopo.h"
#include "CellTopo.h"
#include "CellMesh.h"
#include "BgGrid.h"
#include "GridState.h"
#include "NodeMesh.h"
#include "InterFace.h"
#include "CalcGrid.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

void L2GMapping::CalcL2G( UnsGrid & ggrid, int zid, UnsGrid & grid, G2LMapping & g2l )
{
    this->Alloc( grid );
    this->CalcL2GNode( ggrid, g2l );
    this->CalcL2GFace( ggrid, g2l );
    this->CalcL2GCell( ggrid, zid, g2l );
}

void L2GMapping::Alloc( UnsGrid & grid )
{
    int nNodes = grid.nNodes;
    int nFaces = grid.nFaces;
    int nCells = grid.nCells;

    this->l2g_node.resize( nNodes );
    this->l2g_face.resize( nFaces );
    this->l2g_cell.resize( nCells );
}

void L2GMapping::CalcL2GNode( UnsGrid & ggrid, G2LMapping & g2l )
{
    int nNodes = ggrid.nNodes;

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        if ( g2l.g2l_node[ iNode ] > - 1 )
        {
            this->l2g_node[ g2l.g2l_node[ iNode ] ] = iNode;
        }
    }
}

void L2GMapping::CalcL2GFace( UnsGrid & ggrid, G2LMapping & g2l )
{
    int nFaces = ggrid.nFaces;

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int fid = g2l.g2l_face[ iFace ];
        if ( fid >= 0 )
        {
            this->l2g_face[ fid ] = iFace;
        }
    }
}

void L2GMapping::CalcL2GCell( UnsGrid & ggrid, int zid, G2LMapping & g2l )
{
    int nCells = ggrid.nCells;
    int cid = 0;
    for ( int gcid = 0; gcid < nCells; ++ gcid )
    {
        if ( g2l.gc2lzone[ gcid ] == zid )
        {
            this->l2g_cell[ cid ] = gcid;
            ++ cid;
        }
    }
}

G2LMapping::G2LMapping( UnsGrid & ggrid, int npartprocIn ) : npartproc( npartprocIn )
{
    this->g2l_cell.resize( ggrid.nCells );
    this->g2l_face.resize( ggrid.nFaces );
    this->g2l_node.resize( ggrid.nNodes );
    this->gc2lzone.resize( ggrid.nCells );

}

void G2LMapping::GenerateGC2Z( UnsGrid & ggrid )
{
    if ( npartproc < 2 )
    {
        Fatal( "The number of partitions should be greater than 1!\n" );
    }

    int nCells  = ggrid.nCells;
    int nFaces  = ggrid.nFaces;
    int nBFaces = ggrid.nBFaces;

    std::vector<idx_t> xadj  ( ggrid.nCells + 1 );
    std::vector<idx_t> adjncy( 2 * ( nFaces - nBFaces ) );

    this->GetXadjAdjncy( ggrid, xadj, adjncy );
    this->PartByMetis( nCells, xadj, adjncy );
}
#ifdef ENABLE_METIS
void G2LMapping::GetXadjAdjncy( UnsGrid & ggrid, std::vector<idx_t> & xadj, std::vector<idx_t>& adjncy )
{   
    int  nCells = ggrid.nCells;
    CalcC2C( ggrid );
    LinkField & c2c = ggrid.GetCellMesh().GetCellTopo().c2c;
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

Partition::Partition( const GridConfig & config )
    : sourceFile( config.sourceFile ), partitionType( config.partitionType ), npartproc( config.partitionCount )
{
}

Partition::~Partition()
{
}

void Partition::Run()
{
    UnsGrid & ggrid = this->ReadGrid( this->sourceFile );
    this->GenerateMultiZoneGrid( ggrid );

    ONEFLOW::GenerateMultiZoneCalcGrids( std::move( grids ) );
}

UnsGrid & Partition::ReadGrid( const std::string & sourceFile )
{
    StringField gridFileList;
    gridFileList.push_back( sourceFile );

    Zone::ReadGrid( gridFileList );

    int nZones = ZoneState::nZones;

    if ( nZones > 1 )
    {
        Fatal( " At present, there is no support for multiple blocks such as nZones > 1 !\n" );
    }

    Grid & grid = Zone::GetGridReference();
    UnsGrid * unsGrid = UnsGridCast( &grid );
    if ( ! unsGrid )
    {
        Fatal( "Partition requires an unstructured grid!\n" );
    }
    return *unsGrid;
}

void Partition::GenerateMultiZoneGrid( UnsGrid & ggrid )
{
    this->CreatePart( ggrid );

    this->AllocPart();

    this->BuildCalculationalGrid( ggrid );
}

void Partition::CreatePart( UnsGrid & ggrid )
{
    g2l.emplace( ggrid, this->npartproc );
    g2l->GenerateGC2Z( ggrid );
    this->CalcG2lCell( ggrid );

}

void Partition::AllocPart()
{
    grids.clear();
    grids.resize( static_cast< std::size_t >( npartproc ) );
    for ( int pid = 0; pid < npartproc; ++ pid )
    {
        const int gridType = ONEFLOW::UMESH;
        auto owned = ONEFLOW::CreateGridUnique( gridType );
        owned->level = 0;
        owned->id = pid;
        owned->localId = pid;
        owned->type = gridType;
        grids[ static_cast< std::size_t >( pid ) ] = std::move( owned );
    }
}

void Partition::BuildCalculationalGrid( UnsGrid & ggrid )
{
    for ( int pid = 0; pid < npartproc; ++ pid )
    {
        std::cout << "BuildCalculationalGrid pid = " << pid << " npartproc = " << npartproc << "\n";
        //for unstructured grid, each processor only contains one zone, so pid equal to zid
        this->BuildCalculationalGrid( ggrid, pid );
    }
}

void Partition::CalcG2lCell( UnsGrid & ggrid )
{
    UnsGrid & grid = ggrid;
    G2LMapping & mapping = this->g2l.value();

    IntField zCount( npartproc, 0 );

    int nCells = grid.nCells;
    for ( int cid = 0; cid < nCells; ++ cid )
    {
        int zid = mapping.gc2lzone[ cid ];
        mapping.g2l_cell[ cid ] = zCount[ zid ] ++;
    }
}

void Partition::BuildCalculationalGrid( UnsGrid & ggrid, int zid )
{
    UnsGrid & grid = static_cast< UnsGrid & >( GridAt( grids, zid ) );

    grid.nCells = this->GetNCell( ggrid, zid );

    this->CalcG2lFace( ggrid, zid, grid );
    this->CalcG2lNode( ggrid, grid );

    this->CreateL2g( ggrid, zid, grid );
    this->SetCoor  ( ggrid, grid );
    this->SetGeometricRelationship( ggrid, zid, grid );
}

void Partition::CalcG2lFace( UnsGrid & ggrid, int zid, UnsGrid & grid )
{
    G2LMapping & mapping = this->g2l.value();

    int nCells  = ggrid.nCells;
    int nFaces  = ggrid.nFaces;
    int nBFaces = ggrid.nBFaces;

    IntField & glCell = ggrid.GetFaceTopo().GetLeftCells();
    IntField & grCell = ggrid.GetFaceTopo().GetRightCells();

    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        mapping.g2l_face[ fid ] = - 2;
    }

    //set all face in iZone to -1
    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        int glc = glCell[ fid ];
        int grc = grCell[ fid ];
        if ( mapping.gc2lzone[ glc ] == zid )
        {
            mapping.g2l_face[ fid ] = - 1;
        }
        else if ( grc < nCells && mapping.gc2lzone[ grc ] == zid )
        {
            mapping.g2l_face[ fid ] = - 1;
        }
    }

    int nFaceNow = 0;

    //physical boundary
    for ( int fid = 0; fid < nBFaces; ++ fid )
    {
        if ( mapping.g2l_face[ fid ] == - 1 )
        {
            mapping.g2l_face[ fid ] = nFaceNow ++;
        }
    }

    int nIFaceNow = 0;
    for ( int fid = nBFaces; fid < nFaces; ++ fid )
    {
        if ( mapping.g2l_face[ fid ] == - 1 )
        {
            int glc = glCell[ fid ];
            int grc = grCell[ fid ];
            if ( mapping.gc2lzone[ glc ] != mapping.gc2lzone[ grc ] )
            {
                mapping.g2l_face[ fid ] = nFaceNow ++;
                nIFaceNow ++;
            }
        }
    }

    //inner boundary
    int nBFaceNow = nFaceNow;
    for ( int fid = nBFaces; fid < nFaces; ++ fid )
    {
        if ( mapping.g2l_face[ fid ] == - 1 )
        {
            int glc = glCell[ fid ];
            int grc = grCell[ fid ];

            if ( mapping.gc2lzone[ glc ] == mapping.gc2lzone[ grc ] )
            {
                mapping.g2l_face[ fid ] = nFaceNow ++;
            }
        }
    }

    grid.nFaces  = nFaceNow;
    grid.nBFaces = nBFaceNow;

    InterFace & interFace = *grid.interFace;
    interFace.Set( nIFaceNow );
    grid.nIFaces = nIFaceNow;
}

void Partition::CalcG2lNode( UnsGrid & ggrid, UnsGrid & grid )
{
    G2LMapping & mapping = this->g2l.value();

    int nFaces = ggrid.nFaces;
    int nNodes = ggrid.nNodes;

    LinkField & f2n = ggrid.GetFaceTopo().GetFaces();

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        mapping.g2l_node[ iNode ] = - 2;
    }

    //set iZone g2l->g2l_node to -1
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( mapping.g2l_face[ iFace ] > - 1 )
        {
            int nFNode = f2n[ iFace ].size();
            for ( int iNode = 0; iNode < nFNode; ++ iNode )
            {
                mapping.g2l_node[ f2n[ iFace ][ iNode ] ] = - 1;
            }
        }
    }

    int nLNode = 0;
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        if ( mapping.g2l_node[ iNode ] == - 1 )
        {
            mapping.g2l_node[ iNode ] = nLNode ++;
        }
    }

    grid.nNodes = nLNode;
}

int Partition::GetNCell( UnsGrid & ggrid, int zid )
{
    G2LMapping & mapping = this->g2l.value();

    int nCells = ggrid.nCells;
    int iCount = 0;
    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        if ( mapping.gc2lzone[ iCell ] == zid )
        {
            iCount ++;
        }
    }
    return iCount;
}

void Partition::CreateL2g( UnsGrid & ggrid, int zid, UnsGrid & grid )
{
    this->l2g.CalcL2G( ggrid, zid, grid, this->g2l.value() );
}

void Partition::SetCoor( UnsGrid & ggrid, UnsGrid & grid )
{
    G2LMapping & mapping = this->g2l.value();

    int nNodes = grid.nNodes;
    grid.nodeMesh->CreateNodes( nNodes );

    int iCount = 0;
    for ( int iNode = 0; iNode < ggrid.nNodes; ++ iNode )
    {
        if ( mapping.g2l_node[ iNode ] > - 1 )
        {
            grid.nodeMesh->xN[ iCount ] = ggrid.nodeMesh->xN[ iNode ];
            grid.nodeMesh->yN[ iCount ] = ggrid.nodeMesh->yN[ iNode ];
            grid.nodeMesh->zN[ iCount ] = ggrid.nodeMesh->zN[ iNode ];
            ++ iCount;
        }
    }

    if ( iCount != nNodes )
    {
        std::cout << "error in Partition::SetCoor\n";
    }
}

void Partition::SetGeometricRelationship( UnsGrid & ggrid, int zid, UnsGrid & grid )
{
    this->CalcF2N( ggrid, grid );
    this->SetF2CAndBC( ggrid, zid, grid );
    this->SetInterface( ggrid, zid, grid, this->partitionType );
}

void Partition::CalcF2N( UnsGrid & ggrid, UnsGrid & grid )
{
    G2LMapping & mapping = this->g2l.value();

    LinkField & f2n = grid.GetFaceTopo().GetFaces();
    LinkField & gf2n = ggrid.GetFaceTopo().GetFaces();

    int nFaces = grid.nFaces;
    f2n.resize( nFaces );

    for ( int fid = 0; fid < nFaces; ++ fid )
    {
        int gfid = this->l2g.l2g_face[ fid ];

        int nFNode = gf2n[ gfid ].size();

        for ( int iNode = 0; iNode < nFNode; ++ iNode )
        {
            int gnid = gf2n[ gfid ][ iNode ];
            int nid  = mapping.g2l_node[ gnid ];
            f2n[ fid ].push_back( nid );
        }
    }
}

void Partition::SetF2CAndBC( UnsGrid & ggrid, int zid, UnsGrid & grid )
{
    G2LMapping & mapping = this->g2l.value();

    int nGBFace = ggrid.nBFaces;

    IntField & glCell = ggrid.GetFaceTopo().GetLeftCells();
    IntField & grCell = ggrid.GetFaceTopo().GetRightCells();

    IntField & gbcType = ggrid.GetFaceTopo().GetBcRecord().bcType;

    int nFaces  = grid.nFaces;
    int nBFaces = grid.nBFaces;

    IntField & lCell = grid.GetFaceTopo().GetLeftCells();
    IntField & rCell = grid.GetFaceTopo().GetRightCells();
    lCell.resize( nFaces );
    rCell.resize( nFaces );

    grid.GetFaceTopo().SetNBFaces( nBFaces );

    IntField & local_bcType = grid.GetFaceTopo().GetBcRecord().bcType;

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int gfid = this->l2g.l2g_face[ iFace ];

        int glc = glCell[ gfid ];
        int grc = grCell[ gfid ];

        int lc, rc, bctype;

        if ( gfid < nGBFace )
        {
            rc = - 1;
            lc = mapping.g2l_cell[ glc ];
            bctype = gbcType[ gfid ];
         }
        else
        {
            bctype = -1;
            // int face
            if ( mapping.gc2lzone[ glc ] == zid )
            {
                lc = mapping.g2l_cell[ glc ];
                rc = - 1;
            }
            else if ( mapping.gc2lzone[ grc ] == zid )
            {
                rc = mapping.g2l_cell[ grc ];
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
        int gfid = this->l2g.l2g_face[ iFace ];
        int glc = glCell[ gfid ];
        int grc = grCell[ gfid ];

        int lc = mapping.g2l_cell[ glc ];
        int rc = mapping.g2l_cell[ grc ];

        lCell[ iFace ] = lc;
        rCell[ iFace ] = rc;
    }
}

void Partition::SetInterface( UnsGrid & ggrid, int zid, UnsGrid & grid, int partitionType )
{
    if ( partitionType != 1 ) return;

    G2LMapping & mapping = this->g2l.value();

    InterFace & interFace = *grid.interFace;
    int nIFaces = interFace.nIFaces;
    int nBFaces = grid.nBFaces;

    int nGFace = ggrid.nFaces;
    int nGBFace = ggrid.nBFaces;

    IntField & glCell = ggrid.GetFaceTopo().GetLeftCells();
    IntField & grCell = ggrid.GetFaceTopo().GetRightCells();

    //number of physical boundary face
    int nPBFace = nBFaces - nIFaces;

    for ( int gfid = nGBFace; gfid < nGFace; ++ gfid )
    {
        int fid = mapping.g2l_face[ gfid ];
        if ( fid < nBFaces && fid > - 1 )
        {
            //local interface id
            int ifid = fid - nPBFace;

            int glc = glCell[ gfid ];
            int grc = grCell[ gfid ];

            int leftZone  = mapping.gc2lzone[ glc ];
            int rightZone = mapping.gc2lzone[ grc ];

            int gcid = -1;

            if ( leftZone == zid )
            {
                interFace.idir[ ifid ] = 1;
                gcid = grc;
            }
            else if ( rightZone == zid )
            {
                interFace.idir[ ifid ] = - 1;
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
            interFace.zoneId          [ ifid ] = mapping.gc2lzone[ gcid ];
            interFace.localInterfaceId[ ifid ] = mapping.g2l_cell[ gcid ];
            interFace.localCellId     [ ifid ] = mapping.g2l_cell[ gcid ];
            interFace.i2b             [ ifid ] = mapping.g2l_face[ gfid ];
        }
    }
}

bool FindMatch( UnsGrid & grid, FacePair & facePair )
{
    bool found = false;

    InterFace * interFace = grid.interFace.get();

    if ( ! ONEFLOW::IsValid( interFace ) ) return found;

    int nBFaces = grid.nBFaces;
    int nIFaces = interFace->nIFaces;
    int nPBFace = nBFaces - nIFaces;

    IntField & lCell = grid.GetFaceTopo().GetLeftCells();
    IntField & rCell = grid.GetFaceTopo().GetRightCells();

    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int lc = lCell[ iFace + nPBFace ];
        int rc = rCell[ iFace + nPBFace ];
        int cell_id = MAX( lc, rc );

        if ( ( interFace->zoneId[ iFace ]      == facePair.lf.zone_id ) &&
             ( interFace->localCellId[ iFace ] == facePair.lf.cell_id ) && 
             ( cell_id                         == facePair.rf.cell_id ) )
        {
            interFace->localInterfaceId[ iFace ] = facePair.lf.face_id;
            facePair.rf.face_id = iFace;
            found = true;
            break;
        } 
    }

    return found;
}

EndNameSpace
