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
#include "Grid.h"
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include "Cmpi.h"
#include "Geom.h"
#ifdef PRJ_ENABLE_CGNS
#include <cgnslib.h>
#endif

Grid::Grid()
{
    ;
}

Grid::~Grid()
{
    ;
}

BcTypeMap::BcTypeMap()
{
    ;
}

BcTypeMap::~BcTypeMap()
{
    ;
}

void BcTypeMap::Init()
{
    typedef std::pair< int, std::string > IntStringPair;

    this->bc_map.insert( IntStringPair( BCInterface, "BCInterface" ) );
    this->bc_map.insert( IntStringPair( BCInflow, "BCInflow" ) );
    this->bc_map.insert( IntStringPair( BCOutflow, "BCOutflow" ) );
}

std::string BcTypeMap::GetBcName( int bcId )
{
    return this->bc_map[ bcId ];
}

Interface_t::Interface_t()
{
    this->zone = -1;
    this->cell = -1;
    this->ghost_cell = -1;
    this->face = -1;
}

Interface_t::~Interface_t()
{
    ;
}

IFaceBasic::IFaceBasic()
{
    this->global_face_id = -1;
    this->count = 0;
    this->bctype = -1;
}

IFaceBasic::~IFaceBasic()
{
}

int IFaceBasic::GetNeighborZone()  const
{
    return this->face_pair[ this->nid ].zone;
}

int IFaceBasic::GetNeighborFace()  const
{
    return this->face_pair[ this->nid ].face;
}

Interface_t * IFaceBasic::GetCurrentInterface_t()
{
    return & this->face_pair[ 1 - this->nid ];
}

Interface_t * IFaceBasic::GetNeighborInterface_t()
{
    return & this->face_pair[ this->nid ];
}


void FreeMemory( std::map<int, IFaceBasic*>& global_ptmap )
{
    for ( auto it = global_ptmap.begin(); it != global_ptmap.end(); ++ it )
    {
        delete it->second;
    }
    global_ptmap.clear();
}


PhysicalBoundarySolver::PhysicalBoundarySolver()
{
}

PhysicalBoundarySolver::~PhysicalBoundarySolver()
{
}

void PhysicalBoundarySolver::Init()
{
    for ( int i = 0; i < this->physical_bc_list.size(); ++ i )
    {
        IFaceBasic *f = this->physical_bc_list[ i ];
        Interface_t * fp = & f->face_pair[ 0 ];
        this->physical_faceids.push_back( fp->face );
        this->physical_ghost_cells.push_back( fp->ghost_cell );
        this->physical_bctypes.push_back( f->bctype );
    }
}

BoundarySolver::BoundarySolver()
{
    this->interfaceSolver = new InterfaceSolver();
    this->physicalSolver = new PhysicalBoundarySolver();
}

BoundarySolver::~BoundarySolver()
{
    FreeMemory( this->global_face_map );
    delete this->interfaceSolver;
    delete this->physicalSolver;
}

int BoundarySolver::GetNIFace()
{ 
    return interfaceSolver->GetNIFace();
}

void Insert( std::map<int, int> &local_face_map, int local_pt, int global_pt )
{
    auto it = local_face_map.find( local_pt );
    if ( it == local_face_map.end() )
    {
        local_face_map.insert( std::pair<int,int>(local_pt, global_pt) );
    }
}

void Insert( std::map<int, IFaceBasic *> & global_face_map, int global_face_id, int local_face_id, int zone_id, int cell_id, int ghost_cell_id )
{
    auto it = global_face_map.find( global_face_id );
    if ( it == global_face_map.end() )
    {
        IFaceBasic* f = new IFaceBasic();
        f->global_face_id = global_face_id;
        int idx = 0;
        f->face_pair[ idx ].zone = zone_id;
        f->face_pair[ idx ].cell = cell_id;
        f->face_pair[ idx ].ghost_cell = ghost_cell_id;
        f->face_pair[ idx ].face = local_face_id;
        f->count ++;
        global_face_map.insert( std::pair<int,IFaceBasic *>(global_face_id, f) );
    }
    else
    {
        IFaceBasic * f = it->second;
        int idx = 1;
        f->face_pair[ idx ].zone = zone_id;
        f->face_pair[ idx ].cell = cell_id;
        f->face_pair[ idx ].ghost_cell = ghost_cell_id;
        f->face_pair[ idx ].face = local_face_id;
        f->count ++;
    }
}

IFaceBasic * BoundarySolver::FindIFaceBasic( int id )
{
    auto it = global_face_map.find( id );
    if ( it != global_face_map.end() )
    {
        return it->second;
    }
    return 0;
}

void BoundarySolver::InsertVector(std::vector<int>& a, std::vector<int>& b )
{
    a.insert( a.end(), b.begin(), b.end() );
}

void BoundarySolver::FillBCPoints()
{
    this->interfaceSolver->zone_id = this->zoneId;
    for ( auto it = this->local_face_map.begin(); it != this->local_face_map.end(); ++ it )
    {
        int local_face_id = it->first;
        int global_face_id = it->second;
        IFaceBasic * f = FindIFaceBasic( global_face_id );

        if ( f->bctype == BCInterface )
        {
            this->interfaceSolver->interlist.push_back( f );
        }
        else
        {
            this->physicalSolver->physical_bc_list.push_back( f );
        }
    }
    this->interfaceSolver->Init();
    this->physicalSolver->Init();

    BcTypeMap bcTypeMap;
    bcTypeMap.Init();

    InsertVector( this->bctypes, this->physicalSolver->physical_bctypes );
    InsertVector( this->bctypes, this->interfaceSolver->interface_bctypes );
    InsertVector( this->bc_faceids, this->physicalSolver->physical_faceids );
    InsertVector( this->bc_faceids, this->interfaceSolver->interface_faceids );
    InsertVector( this->bc_ghostcells, this->physicalSolver->physical_ghost_cells );
    InsertVector( this->bc_ghostcells, this->interfaceSolver->interface_ghost_cells );

    //std::printf( "bcinfo......\n" );
    //std::printf("local_face_map.size()=%zd\n", local_face_map.size() );
    //std::printf("local_face_map.size()=%zd\n", local_face_map.size() );
    for ( int iface = 0; iface < this->bc_faceids.size(); ++ iface )
    {
        int bc_faceid = this->bc_faceids[ iface ];
        int bctype = this->bctypes[ iface ];
        std::string bcName = bcTypeMap.GetBcName( bctype );
        //std::printf("bc_faceid=%d bctype=%d, bcName=%s\n", bc_faceid, bctype, bcName.c_str() );
    }
    //std::printf("\n");

}

void BoundarySolver::PrintBcInfo()
{
    std::printf( "BoundarySolver::Init pt: zoneId = %d nZones = %d\n", zoneId, nZones );
    std::printf("global face map : \n" );
    for ( auto it = global_face_map.begin(); it != global_face_map.end(); ++ it )
    {
        IFaceBasic *f = it->second;
        Interface_t * f1 = &f->face_pair[ 0 ];
        Interface_t * f2 = &f->face_pair[ 1 ];
        std::printf( "%d->%d ", it->first, f->count );
        std::printf( "(zone,cell,ghost)[(%d,%d,%d),(%d,%d,%d)] bcType = %d\n", \
            f1->zone, f1->cell, f1->ghost_cell, \
            f2->zone, f2->cell, f2->ghost_cell, \
            f->bctype );
    }
    std::printf("\n");
    std::printf("local face map : \n" );
    for ( auto it = local_face_map.begin(); it != local_face_map.end(); ++ it )
    {
        std::printf( "%d -> %d \n", it->first, it->second );
    }
    std::printf("\n");
}

void BoundarySolver::MarkInterface()
{
    for ( auto it = global_face_map.begin(); it != global_face_map.end(); ++ it )
    {
        if ( it->second->count == 2 )
        {
            it->second->bctype = BCInterface;
        }
    }
}

void BoundarySolver::SetBcType( int bc_face_id, int bcType )
{
    auto it = global_face_map.find( bc_face_id );
    if ( it != global_face_map.end() )
    {
        it->second->bctype = bcType;
    }
}

void BoundarySolver::Init( int zoneId, int nZones, int ni )
{
    this->zoneId = zoneId;
    this->nZones = nZones;
    this->ni = ni;
    //local pt set : 30 40
    //global pt map : 0->1 10->2 20->2 30->2 40->1
    int ishift = 0;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        int ni_local = Geom_t::zone_nis[ iZone ];
        int local_face_id0 = 1 + 0;
        int local_face_id1 = 1 + ni_local - 1;
        int global_face_id0 = ishift + 1 + 0;
        int global_face_id1 = ishift + 1 + ni_local - 1;
        ishift += ( ni_local - 1 );
        int cell_id0 = local_face_id0 + 1;
        int cell_id1 = local_face_id1 - 1;
        int ghost_cell_id0 = local_face_id0 - 1;
        int ghost_cell_id1 = local_face_id1 + 1;
        Insert( global_face_map, global_face_id0, local_face_id0, iZone, cell_id0, ghost_cell_id0 );
        Insert( global_face_map, global_face_id1, local_face_id1, iZone, cell_id1, ghost_cell_id1 );

        if ( iZone == zoneId )
        {
            Insert( local_face_map, local_face_id0, global_face_id0 );
            Insert( local_face_map, local_face_id1, global_face_id1 );
        }
    }
    this->MarkInterface();

    int left_bcface_id = 1;
    int right_bcface_id = 1 + ishift;
    std::printf( " right_bcface_id =%d \n", right_bcface_id );
    this->SetBcType( left_bcface_id, BCInflow );
    this->SetBcType( right_bcface_id, BCOutflow );

    this->FillBCPoints();
}


IData::IData()
{
    ;
}

IData::~IData()
{
    ;
}

void IData::Print()
{
    std::printf( " IData::Print() zone_id =%d \n", zone_id );
    for ( int i = 0; i < this->interface_neighbor_faceids.size(); ++ i )
    {
        std::printf( "%d ", this->interface_neighbor_faceids[ i ] );
    }
    std::printf( "\n" );
}

void IData::SetSendQs( float * q )
{
    int nSize = this->send_cells.size();
    this->send_qs.resize( nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        int id = this->send_cells[i];
        send_qs[ i ] = q[ id ];
    }
}

void IData::GetRecvQs( float * q )
{
    int nSize = this->recv_cells.size();
    this->recv_qs.resize( nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        int id = this->send_cells[i];
        recv_qs[ i ] = q[ id ];
    }
}

InterfaceSolver::InterfaceSolver()
{
    ;
}

InterfaceSolver::~InterfaceSolver()
{
    ;
}

bool CmpIFaceBasicPtr(const IFaceBasic * a, const IFaceBasic * b)
{
    int ia = 1 - a->nid;
    int ib = 1 - b->nid;
    return a->face_pair[ ia ].face < b->face_pair[ ib ].face;
}

bool CmpIFaceBasicNeighborPtr(const IFaceBasic * a, const IFaceBasic * b)
{
    int a_neighbor_zone = a->GetNeighborZone();
    int b_neighbor_zone = b->GetNeighborZone();
    if ( a_neighbor_zone != b_neighbor_zone )
    {
        return a_neighbor_zone < b_neighbor_zone;
    }
    return a->GetNeighborFace() < b->GetNeighborFace();
}


void InterfaceSolver::SortInterface()
{
    //std::printf( " SetNeighborData nSize =%zd \n", this->interlist.size() );
    for ( int i = 0; i < this->interlist.size(); ++ i )
    {
        IFaceBasic *f = this->interlist[ i ];
        if ( f->face_pair[0].zone == this->zone_id )
        {
            f->nid = 1;
        }
        else
        {
            f->nid = 0;
        }
    }
    std::sort( this->interlist.begin(), this->interlist.end(), CmpIFaceBasicNeighborPtr );
    //for ( int i = 0; i < this->interlist.size(); ++ i )
    //{
    //    IFaceBasic *f = this->interlist[ i ];
    //    std::printf( " i =%d neighbor_zone = %d neighbor_face = %d \n", i, f->face_pair[f->nid].zone, f->face_pair[f->nid].face);
    //}
}

void InterfaceSolver::SetInterface()
{
    for ( int i = 0; i < this->interlist.size(); ++ i )
    {
        IFaceBasic * f = this->interlist[ i ];
        Interface_t * finter = & f->face_pair[ 1 - f->nid ];
        Interface_t * fneighbor = & f->face_pair[ f->nid ];
        this->interface_bctypes.push_back( f->bctype );
        this->interface_faceids.push_back( finter->face );
        this->interface_ghost_cells.push_back( finter->ghost_cell );
        this->interface_neighbor_zoneids.push_back( fneighbor->zone );
        this->interface_neighbor_faceids.push_back( fneighbor->face );
        this->interface_neighbor_ghost_cells.push_back( fneighbor->ghost_cell );
        
    }
}

void InterfaceSolver::SetNeighborData()
{
    int nSize = this->interlist.size();
    //std::printf( " SetNeighborData nSize =%d \n", nSize );
    if ( nSize == 0 ) return;

    int ist = 0;
    int icount = 0;
    while ( true )
    {
        std::vector<int> interfaces;
        IData iData;
        int neighbor_zone_id_start = this->interlist[ ist ]->GetNeighborZone();
        int ipos = this->neighbor_datas.size();
        iData.zone_id = neighbor_zone_id_start;

        for ( int i = ist; i < this->interlist.size(); ++ i )
        {
            IFaceBasic * f = this->interlist[ i ];
            int neighbor_zone_now = f->GetNeighborZone();
            if ( neighbor_zone_id_start != neighbor_zone_now )
            {
                ist = i;
                break;
            }
            int jpos = iData.interface_neighbor_faceids.size();
            int neighbor_face = f->GetNeighborFace();

            this->interface_ipos.push_back( ipos );
            this->interface_jpos.push_back( jpos );
            iData.interface_neighbor_faceids.push_back( neighbor_face );
            iData.tmplist.push_back( f );
            Interface_t * ft = f->GetCurrentInterface_t();
            iData.recv_cells.push_back( ft->cell );
            icount ++;
        }
        this->neighbor_datas.push_back( iData );
        if ( icount == this->interlist.size() ) break;
    }
    ReorderCurrentData();
}

void InterfaceSolver::ReorderCurrentData()
{
    int nSize = this->neighbor_datas.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        IData &iData = this->neighbor_datas[ i ];
        std::sort( iData.tmplist.begin(), iData.tmplist.end(), CmpIFaceBasicPtr );
        for ( int i = 0; i < iData.tmplist.size(); ++ i )
        {
            IFaceBasic *f = iData.tmplist[ i ];
            Interface_t * fp = f->GetCurrentInterface_t();
            iData.send_cells.push_back( fp->cell );
            //std::printf( " i =%d current_zone = %d current_face = %d \n", i, fp->zone, fp->face);
        }
    }
}

float InterfaceSolver::GetBcValue( int iface )
{
    int ipos = this->interface_ipos[ iface ];
    int jpos = this->interface_jpos[ iface ];
    return this->neighbor_datas[ ipos ].recv_qs[ jpos ];
}

void InterfaceSolver::ShowInfo( int iface )
{
    int ipos = this->interface_ipos[ iface ];
    int jpos = this->interface_jpos[ iface ];
    int neighbor_zone_id = this->interface_neighbor_zoneids[ iface ];
    int neighbor_face_id = this->interface_neighbor_faceids[ iface ];
    std::printf( "ShowInfo ipos = %d, jpos = %d, zoneid = %d, neighbor_zone_id = %d, neighbor_face_id = %d \n",\
        ipos, jpos, this->zone_id, neighbor_zone_id, neighbor_face_id );
}

void InterfaceSolver::Init()
{
    //std::printf( " InterfaceSolver::Init()\n" );
    this->SortInterface();
    this->SetInterface();
    this->SetNeighborData();

    //std::printf( " myzoneId =%d neighbor_datas.size() = %zd --------------------\n", zone_id, neighbor_datas.size() );
    for ( int i = 0; i < this->neighbor_datas.size(); ++ i )
    {
        //this->neighbor_datas[ i ].Print();
    }
    this->SwapData();
}

void InterfaceSolver::SwapData()
{
}

void InterfaceSolver::PrepareSendData( float * q )
{
    for ( int i = 0; i < this->neighbor_datas.size(); ++ i )
    {
        IData & idata = this->neighbor_datas[ i ];
        idata.SetSendQs(q);
    }
}

void InterfaceSolver::SwapData( float * q )
{
    PrepareSendData( q );
#ifdef PRJ_ENABLE_MPI
    for ( int i = 0; i < this->neighbor_datas.size(); ++ i )
    {
        IData & idata = this->neighbor_datas[ i ];
        int neighbor_zoneid = idata.zone_id;
        int ip = neighbor_zoneid;
        int tag = 0;
        if ( neighbor_zoneid > this->zone_id )
        {
            int n = idata.send_qs.size();
            MPI_Send(&n, 1, MPI_INT, ip, tag, MPI_COMM_WORLD);
            MPI_Send(idata.send_qs.data(), n, MPI_FLOAT, ip, tag, MPI_COMM_WORLD);
        }
        else
        {
            int n = -1;
            MPI_Recv(&n, 1, MPI_INT, ip, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            idata.recv_qs.resize( n, -1 );
            MPI_Recv(idata.recv_qs.data(), n, MPI_FLOAT, ip, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for ( int i = 0; i < this->neighbor_datas.size(); ++ i )
    {
        IData & idata = this->neighbor_datas[ i ];
        int neighbor_zoneid = idata.zone_id;
        int ip = neighbor_zoneid;
        int tag = 0;
        if ( neighbor_zoneid > this->zone_id )
        {
            int n = -1;
            MPI_Recv(&n, 1, MPI_INT, ip, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            idata.recv_qs.resize( n, -1 );
            MPI_Recv(idata.recv_qs.data(), n, MPI_FLOAT, ip, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            int n = idata.send_qs.size();
            MPI_Send(&n, 1, MPI_INT, ip, tag, MPI_COMM_WORLD);
            MPI_Send(idata.send_qs.data(), n, MPI_FLOAT, ip, tag, MPI_COMM_WORLD);
        }
    }
#endif
}
